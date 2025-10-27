import argparse
import itertools
import math
import os
import time
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from SB_SupportFunction import read_input
from evaluate_all_sols_check import evaluate_all_sols
from linear_programming_scheduler import LinearProgrammingScheduler

VALUE_HEAVY = 40
DEFAULT_PROGRESS_INTERVAL = 10_000


@dataclass
class ExhaustiveResult:
    permutation: List[int]
    tardiness: float
    list_picker: List[List[int]]
    batches: List[List[int]]
    batch_completion: Dict[int, float]
    order_tardiness: Dict[int, float]


@dataclass
class ExhaustiveMetrics:
    total_permutations: int
    evaluated_permutations: int
    elapsed_seconds: float
    completed: bool


class ExhaustivePermutationOptimizer:
    """Enumerate all permutations of the items and evaluate tardiness."""

    def __init__(
        self,
        name_path_input: str,
        *,
        progress_interval: int = DEFAULT_PROGRESS_INTERVAL,
        enumeration_time_limit: Optional[float] = None,
        scheduler_time_limit: Optional[float] = None,
        scheduler_progress_interval: int = 100_000,
        checkpoint_dir: Optional[str] = None,
        max_items: Optional[int] = None,
    ) -> None:
        self.name_path_input = name_path_input
        self.progress_interval = max(1, progress_interval)
        self.enumeration_time_limit = enumeration_time_limit
        self.checkpoint_dir = checkpoint_dir

        df_item_pool, _ = read_input(name_path_input)
        if max_items is not None and 0 < max_items < len(df_item_pool):
            df_item_pool = df_item_pool.iloc[:max_items]
            df_item_pool = df_item_pool.reset_index()
            if 'index' in df_item_pool.columns:
                df_item_pool = df_item_pool.drop(['index'], axis=1)

        self.df_item_pool = df_item_pool
        self.num_items = self.df_item_pool.shape[0]
        self.indices: List[int] = list(range(self.num_items))

        self.heavy_item_set = set(
            idx for idx, value in enumerate(self.df_item_pool['weight']) if float(value) >= VALUE_HEAVY
        )

        self.scheduler = LinearProgrammingScheduler(
            max_seconds=scheduler_time_limit, progress_interval=scheduler_progress_interval
        )

        self.total_permutations = math.factorial(self.num_items) if self.num_items > 0 else 0
        self._start_time = 0.0
        self._evaluated = 0
        self._timed_out = False

        self.best_result: Optional[ExhaustiveResult] = None

    def _should_stop(self) -> bool:
        if self.enumeration_time_limit is None:
            return False
        elapsed = time.perf_counter() - self._start_time
        return elapsed >= self.enumeration_time_limit

    def _iter_permutations(self) -> Iterable[Sequence[int]]:
        return itertools.permutations(self.indices)

    def _save_checkpoint(self) -> None:
        if self.checkpoint_dir is None or self.best_result is None:
            return
        os.makedirs(self.checkpoint_dir, exist_ok=True)
        checkpoint_path = os.path.join(self.checkpoint_dir, 'exhaustive_best.txt')
        with open(checkpoint_path, 'w', encoding='utf-8') as handle:
            handle.write('Best tardiness: {:.6f}\n'.format(self.best_result.tardiness))
            handle.write('Best permutation: {}\n'.format(' '.join(map(str, self.best_result.permutation))))
            handle.write('Permutations evaluated: {}\n'.format(self._evaluated))
            handle.write('Total permutations: {}\n'.format(self.total_permutations))

    def _report_progress(self) -> None:
        if self._evaluated == 0:
            return
        elapsed = time.perf_counter() - self._start_time
        rate = self._evaluated / elapsed if elapsed > 0 else 0.0
        remaining = self.total_permutations - self._evaluated
        est_remaining = remaining / rate if rate > 0 else float('inf')
        best = self.best_result.tardiness if self.best_result is not None else float('inf')
        print(
            '[Exhaustive] evaluated {}/{} permutations ({:.6%}); elapsed {:.2f}s; {:.2f} perms/s; best tardiness {:.3f}; est. remaining {:.2f}s'.format(
                f'{self._evaluated:,}',
                f'{self.total_permutations:,}',
                self._evaluated / self.total_permutations if self.total_permutations else 0.0,
                elapsed,
                rate,
                best,
                est_remaining,
            )
        )

    def solve(self) -> ExhaustiveResult:
        if self.num_items == 0:
            raise RuntimeError('No items available for enumeration')

        self._start_time = time.perf_counter()
        for permutation in self._iter_permutations():
            if self._should_stop():
                self._timed_out = True
                break

            self._evaluated += 1
            list_picker, total_tardiness, batches = evaluate_all_sols(
                list(permutation),
                self.df_item_pool.copy(),
                self.heavy_item_set,
                self.name_path_input,
                scheduler=self.scheduler,
            )

            tardiness_value = float(total_tardiness)
            if self.best_result is None or tardiness_value < self.best_result.tardiness:
                self.best_result = ExhaustiveResult(
                    permutation=list(permutation),
                    tardiness=tardiness_value,
                    list_picker=[seq[:] for seq in list_picker],
                    batches=[batch[:] for batch in batches],
                    batch_completion=self.scheduler.last_batch_completion.copy(),
                    order_tardiness=self.scheduler.last_order_tardiness.copy(),
                )
                self._save_checkpoint()

            if self._evaluated % self.progress_interval == 0:
                self._report_progress()

        if self._evaluated % self.progress_interval != 0:
            self._report_progress()

        if self.best_result is None:
            raise RuntimeError('Enumeration did not produce any feasible solution')

        return self.best_result

    def metrics(self) -> ExhaustiveMetrics:
        elapsed = time.perf_counter() - self._start_time if self._start_time else 0.0
        return ExhaustiveMetrics(
            total_permutations=self.total_permutations,
            evaluated_permutations=self._evaluated,
            elapsed_seconds=elapsed,
            completed=not self._timed_out and self._evaluated == self.total_permutations,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description='Exhaustive permutation optimizer for tardiness minimisation')
    parser.add_argument('instance', nargs='?', default='1R-20I-150C-2P', help='Dataset folder name')
    parser.add_argument('--progress-interval', type=int, default=DEFAULT_PROGRESS_INTERVAL, help='Permutations between progress reports')
    parser.add_argument('--time-limit', type=float, default=None, help='Optional wall-clock limit for permutation enumeration (seconds)')
    parser.add_argument('--scheduler-time-limit', type=float, default=None, help='Optional wall-clock limit for each scheduling solve')
    parser.add_argument('--scheduler-progress-interval', type=int, default=100_000, help='Nodes between scheduling progress reports')
    parser.add_argument('--checkpoint-dir', type=str, default=None, help='Directory for writing the best-so-far checkpoint')
    parser.add_argument('--max-items', type=int, default=None, help='Limit the enumeration to the first N items (for experimentation/testing)')
    args = parser.parse_args()

    optimizer = ExhaustivePermutationOptimizer(
        args.instance,
        progress_interval=args.progress_interval,
        enumeration_time_limit=args.time_limit,
        scheduler_time_limit=args.scheduler_time_limit,
        scheduler_progress_interval=args.scheduler_progress_interval,
        checkpoint_dir=args.checkpoint_dir,
        max_items=args.max_items,
    )

    print('Exhaustive permutation optimizer')
    print('--------------------------------')
    print(f'Instance: {optimizer.name_path_input}')
    print(f'Items: {optimizer.num_items}')
    print(f'Total permutations: {optimizer.total_permutations:,}')

    result = optimizer.solve()
    metrics = optimizer.metrics()

    print('\nBest solution summary')
    print('---------------------')
    print(f'Best tardiness: {result.tardiness:.3f}')
    print(f'Permutation: {result.permutation}')
    print('Picker schedules:')
    for idx, sequence in enumerate(result.list_picker, start=1):
        formatted = ', '.join(str(batch) for batch in sequence)
        print(f'  Picker {idx}: [{formatted}]')

    print('\nBatch completion times:')
    for batch_id in sorted(result.batch_completion):
        completion = result.batch_completion[batch_id]
        print(f'  Batch {batch_id}: {completion:.3f}')

    print('\nOrder tardiness:')
    for order_id in sorted(result.order_tardiness):
        tardiness = result.order_tardiness[order_id]
        print(f'  Order {order_id}: {tardiness:.3f}')

    print('\nEnumeration metrics')
    print('-------------------')
    print(f'Permutations evaluated: {metrics.evaluated_permutations:,} / {metrics.total_permutations:,}')
    print(f'Elapsed time: {metrics.elapsed_seconds:.2f} seconds')
    status = 'complete' if metrics.completed else 'partial (time limit reached)'
    print(f'Status: {status}')


if __name__ == '__main__':
    main()
