from dataclasses import dataclass
from typing import Dict, List

from SB_SupportFunction import read_input
from evaluate_all_sols_check import evaluate_all_sols
from linear_programming_scheduler import LinearProgrammingScheduler

VALUE_HEAVY = 40
RANDOM_SEED = 2222


@dataclass
class LinearProgrammingResult:
    permutation: List[int]
    tardiness: float
    list_picker: List[List[int]]
    batches: List[List[int]]
    batch_completion: Dict[int, float]
    order_tardiness: Dict[int, float]


class LinearProgrammingOptimizer:
    def __init__(self, name_path_input: str):
        self.name_path_input = name_path_input
        self.df_item_pool, _ = read_input(name_path_input)
        self.num_items = self.df_item_pool.shape[0]
        self.columns = list(self.df_item_pool.columns)
        self.heavy_item_set = set(self.df_item_pool[self.df_item_pool['weight'] >= VALUE_HEAVY].index)
        self.scheduler = LinearProgrammingScheduler()
        self.best_result: LinearProgrammingResult | None = None

    def _value(self, index: int, column: str) -> float:
        col_idx = self.columns.index(column)
        return self.df_item_pool.iloc[index][col_idx]

    def _candidate_permutations(self) -> List[List[int]]:
        indices = list(range(self.num_items))
        candidates: List[List[int]] = []

        candidates.append(indices)

        due_sorted = sorted(indices, key=lambda idx: (self._value(idx, 'duedate'), idx))
        candidates.append(due_sorted)
        candidates.append(list(reversed(due_sorted)))

        weight_sorted = sorted(indices, key=lambda idx: (-self._value(idx, 'weight'), idx))
        candidates.append(weight_sorted)

        order_sorted = sorted(indices, key=lambda idx: (self._value(idx, 'order'), self._value(idx, 'duedate'), idx))
        candidates.append(order_sorted)

        heavy_first = sorted(indices, key=lambda idx: (0 if idx in self.heavy_item_set else 1, self._value(idx, 'duedate'), idx))
        candidates.append(heavy_first)

        unique: List[List[int]] = []
        seen = set()
        for candidate in candidates:
            key = tuple(candidate)
            if key in seen:
                continue
            seen.add(key)
            unique.append(candidate)
        return unique

    def solve(self) -> LinearProgrammingResult:
        best_tardiness = float('inf')
        best_result: LinearProgrammingResult | None = None

        for permutation in self._candidate_permutations():
            list_picker, total_tardiness, batches = evaluate_all_sols(
                permutation,
                self.df_item_pool.copy(),
                self.heavy_item_set,
                self.name_path_input,
                scheduler=self.scheduler,
            )

            tardiness_value = float(total_tardiness)
            if tardiness_value < best_tardiness:
                best_tardiness = tardiness_value
                best_result = LinearProgrammingResult(
                    permutation=permutation.copy(),
                    tardiness=tardiness_value,
                    list_picker=[seq[:] for seq in list_picker],
                    batches=[batch[:] for batch in batches],
                    batch_completion=self.scheduler.last_batch_completion.copy(),
                    order_tardiness=self.scheduler.last_order_tardiness.copy(),
                )

        if best_result is None:
            raise RuntimeError("Failed to compute a feasible solution")

        self.best_result = best_result
        return best_result


def main() -> None:
    optimizer = LinearProgrammingOptimizer('1R-20I-150C-2P')
    result = optimizer.solve()

    print('Linear Programming Optimizer')
    print('----------------------------')
    print(f'Instance: {optimizer.name_path_input}')
    print(f'Best tardiness: {result.tardiness:.3f}')
    print(f'Permutation length: {len(result.permutation)} items')
    print('Picker schedules:')
    for idx, sequence in enumerate(result.list_picker, start=1):
        formatted = ', '.join(str(batch) for batch in sequence)
        print(f'  Picker {idx}: [{formatted}]')

    print('Batch completion times:')
    for batch_id in sorted(result.batch_completion):
        completion = result.batch_completion[batch_id]
        print(f'  Batch {batch_id}: {completion:.3f}')

    print('Order tardiness:')
    for order_id in sorted(result.order_tardiness):
        tardiness = result.order_tardiness[order_id]
        print(f'  Order {order_id}: {tardiness:.3f}')
if __name__ == '__main__':
    main()
