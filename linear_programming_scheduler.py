import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Set

import pandas as pd


@dataclass
class ScheduleResult:
    picker_sequences: List[List[int]]
    batch_completion: Dict[int, float]
    order_tardiness: Dict[int, float]
    total_tardiness: float


class _BranchAndBoundSolver:
    def __init__(
        self,
        batch_ids: List[int],
        process_times: Dict[int, float],
        orders_by_batch: Dict[int, Set[int]],
        batches_by_order: Dict[int, Set[int]],
        due_dates: Dict[int, float],
        num_picker: int,
    ) -> None:
        self.batch_ids = batch_ids
        self.process_times = process_times
        self.orders_by_batch = orders_by_batch
        self.batches_by_order = batches_by_order
        self.due_dates = due_dates
        self.num_picker = num_picker

        self.best_tardiness: Optional[float] = None
        self.best_sequences: Optional[List[List[int]]] = None

    def solve(self) -> ScheduleResult:
        machine_times = [0.0 for _ in range(self.num_picker)]
        sequences: List[List[int]] = [[] for _ in range(self.num_picker)]
        current_completion: Dict[int, float] = {}
        pending_orders: Set[int] = set(self.due_dates.keys())

        self._search(
            remaining=set(self.batch_ids),
            machine_times=machine_times,
            sequences=sequences,
            current_completion=current_completion,
            pending_orders=pending_orders,
            partial_tardiness=0.0,
        )

        if self.best_sequences is None:
            raise RuntimeError("Failed to produce a feasible schedule")

        batch_completion: Dict[int, float] = {}
        for seq in self.best_sequences:
            time = 0.0
            for batch in seq:
                time += self.process_times[batch]
                batch_completion[batch] = time

        order_tardiness: Dict[int, float] = {}
        total_tardiness = 0.0
        for order, due_date in self.due_dates.items():
            related_batches = self.batches_by_order.get(order, set())
            if not related_batches:
                order_tardiness[order] = 0.0
                continue
            completion = max(batch_completion[batch] for batch in related_batches)
            tardiness = max(0.0, completion - due_date)
            order_tardiness[order] = tardiness
            total_tardiness += tardiness

        return ScheduleResult(
            picker_sequences=[seq[:] for seq in self.best_sequences],
            batch_completion=batch_completion,
            order_tardiness=order_tardiness,
            total_tardiness=total_tardiness,
        )

    def _search(
        self,
        *,
        remaining: Set[int],
        machine_times: List[float],
        sequences: List[List[int]],
        current_completion: Dict[int, float],
        pending_orders: Set[int],
        partial_tardiness: float,
    ) -> None:
        if not remaining:
            if self.best_tardiness is None or partial_tardiness < self.best_tardiness:
                self.best_tardiness = partial_tardiness
                self.best_sequences = [seq[:] for seq in sequences]
            return

        if self.best_tardiness is not None and partial_tardiness >= self.best_tardiness:
            return

        for batch in sorted(remaining):
            orders = self.orders_by_batch.get(batch, set())
            for picker_idx in range(self.num_picker):
                if picker_idx > 0 and not sequences[picker_idx] and all(
                    len(sequences[prev]) == 0 for prev in range(picker_idx)
                ):
                    continue

                start_time = machine_times[picker_idx]
                completion_time = start_time + self.process_times[batch]

                new_machine_times = machine_times[:]
                new_machine_times[picker_idx] = completion_time

                new_sequences = [seq[:] for seq in sequences]
                new_sequences[picker_idx].append(batch)

                new_current_completion = current_completion.copy()
                new_pending_orders = pending_orders.copy()
                new_partial = partial_tardiness

                for order in orders:
                    prev_completion = new_current_completion.get(order, 0.0)
                    if completion_time > prev_completion:
                        new_current_completion[order] = completion_time

                new_remaining = remaining.copy()
                new_remaining.remove(batch)

                completed_orders: List[int] = []
                for order in orders:
                    if order not in new_pending_orders:
                        continue
                    if not (self.batches_by_order.get(order, set()) & new_remaining):
                        completed_orders.append(order)

                for order in completed_orders:
                    new_pending_orders.remove(order)
                    tardiness = max(0.0, new_current_completion.get(order, 0.0) - self.due_dates[order])
                    new_partial += tardiness

                if self.best_tardiness is not None and new_partial >= self.best_tardiness:
                    continue

                self._search(
                    remaining=new_remaining,
                    machine_times=new_machine_times,
                    sequences=new_sequences,
                    current_completion=new_current_completion,
                    pending_orders=new_pending_orders,
                    partial_tardiness=new_partial,
                )


class LinearProgrammingScheduler:
    def __init__(self) -> None:
        self.last_sequences: List[List[int]] = []
        self.last_batch_completion: Dict[int, float] = {}
        self.last_order_tardiness: Dict[int, float] = {}
        self.last_total_tardiness: float = math.inf

    def __call__(self, df_item_pool: pd.DataFrame, num_picker: int):
        rows = df_item_pool.to_dicts()
        batches: Dict[int, Dict[str, object]] = {}
        order_due_dates: Dict[int, float] = {}

        for idx, row in enumerate(rows):
            batch_id = int(row.get('batch', 0))
            if batch_id == 0:
                continue
            order_id = int(row['order'])
            due_date = float(row['duedate'])
            if order_id not in order_due_dates or due_date < order_due_dates[order_id]:
                order_due_dates[order_id] = due_date

            info = batches.setdefault(
                batch_id,
                {'items': [], 'orders': set(), 'process_time': float(row.get('process_time', 0.0))},
            )
            info['items'].append(idx)
            info['orders'].add(order_id)
            if row.get('process_time') is not None:
                info['process_time'] = float(row['process_time'])

        batch_ids = sorted(batches.keys())
        process_times = {batch_id: float(batches[batch_id]['process_time']) for batch_id in batch_ids}
        orders_by_batch = {batch_id: set(batches[batch_id]['orders']) for batch_id in batch_ids}
        batches_by_order: Dict[int, Set[int]] = {}
        for batch_id, orders in orders_by_batch.items():
            for order in orders:
                batches_by_order.setdefault(order, set()).add(batch_id)

        solver = _BranchAndBoundSolver(
            batch_ids=batch_ids,
            process_times=process_times,
            orders_by_batch=orders_by_batch,
            batches_by_order=batches_by_order,
            due_dates=order_due_dates,
            num_picker=num_picker,
        )

        result = solver.solve()

        self.last_sequences = result.picker_sequences
        self.last_batch_completion = result.batch_completion
        self.last_order_tardiness = result.order_tardiness
        self.last_total_tardiness = result.total_tardiness

        df_item_pool['CompletionTime'] = 0.0
        df_item_pool['TardinessOrder'] = 0.0
        for batch_id, completion in result.batch_completion.items():
            df_item_pool.loc[df_item_pool['batch'] == batch_id, 'CompletionTime'] = float(format(completion, '.3f'))
        for order_id, tardiness in result.order_tardiness.items():
            df_item_pool.loc[df_item_pool['order'] == order_id, 'TardinessOrder'] = float(format(tardiness, '.3f'))

        max_order = int(df_item_pool['order'].max()) if hasattr(df_item_pool, 'max') else 0
        list_tardiness_each_order = [0.0 for _ in range(max_order + 1)]
        for order_id, tardiness in result.order_tardiness.items():
            if order_id < len(list_tardiness_each_order):
                list_tardiness_each_order[order_id] = float(format(tardiness, '.3f'))

        return (
            [seq[:] for seq in result.picker_sequences],
            list_tardiness_each_order,
            float(format(result.total_tardiness, '.3f')),
            df_item_pool,
        )
