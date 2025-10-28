import csv
import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from calculate_process_time import calculate_completion_time
from batching_item_check import batching_open
from sequencing_assignment_algorithms import ESDR_algorithms
from routing import combined_routing, precedence_constrained_routing, s_shape_routing


@dataclass
class BatchPreparation:
    dataframe: pd.DataFrame
    batches: List[List[int]]
    index_items_in_batch: List[List[int]]


def _parse_setting_value(value: str) -> Any:
    try:
        if value is None:
            return value
        if '.' in value:
            return float(value)
        return int(value)
    except ValueError:
        return value


def _prepare_batches(new_pop_sol, df_item_poor_batch, heavy_item_set, settings: Dict[str, Any], name_path_input: str) -> BatchPreparation:
    record_list_batch_PSO: List[List[int]] = []
    list_batching_item_PSO: List[List[int]] = []
    list_index_item_in_batch: List[List[int]] = []
    list_index = new_pop_sol

    method_batch = int(settings.get('method_batch', 2))
    capacity_picker_1 = int(settings.get('capacity_picker', 0))
    value_threshold_1 = int(settings.get('value_threshold', 0))

    if method_batch != 2:
        raise ValueError("Only method_batch == 2 is supported in the current implementation")

    list_batching_item_PSO, df_item_poor_new_PSO, list_index_item_in_batch = batching_open(
            df_item_poor_batch,
            list_index,
            capacity_picker_1,
            value_threshold_1,
            heavy_item_set,
            name_path_input,
        )

    def _coerce_nested(values: List[List[Any]]) -> List[List[int]]:
        coerced: List[List[int]] = []
        for batch in values:
            new_batch: List[int] = []
            for item in batch:
                if hasattr(item, 'tolist'):
                    as_list = item.tolist()
                    if isinstance(as_list, list):
                        new_batch.append(int(as_list[0]))
                    else:
                        new_batch.append(int(as_list))
                else:
                    new_batch.append(int(item))
            coerced.append(new_batch)
        return coerced

    list_batching_item_PSO = _coerce_nested(list_batching_item_PSO)
    if list_index_item_in_batch:
        list_index_item_in_batch = _coerce_nested(list_index_item_in_batch)

    list_batching = []
    num_batch = len(list_batching_item_PSO)
    for i in range(num_batch):
        for j in list_batching_item_PSO[i]:
            list_batching.append(j)

    if list_batching_item_PSO:
        list_batching = sum(list_batching_item_PSO, [])
        record_list_batch_PSO.append(list_batching)

    method_routing = int(settings.get('method_routing', 3))
    aisle = int(settings.get('aisle', 0))
    each_aisle_item = int(settings.get('each_aisle_item', 0))
    rack_x = int(settings.get('rack_x', 0))
    rack_y = int(settings.get('rack_y', 0))
    aisle_x = int(settings.get('aisle_x', 0))
    enter_aisle = int(settings.get('enter_aisle', 0))
    distance_y = int(settings.get('distance_y', 0))

    list_distance_batch: List[float] = []
    for item_list in list_batching_item_PSO:
        if method_routing == 1:
            distance_each_batch = s_shape_routing(
                item_list,
                aisle,
                each_aisle_item,
                rack_x,
                rack_y,
                aisle_x,
                enter_aisle,
                distance_y,
            )
        elif method_routing == 2:
            distance_each_batch = combined_routing(
                item_list,
                aisle,
                each_aisle_item,
                rack_x,
                rack_y,
                aisle_x,
                enter_aisle,
                distance_y,
            )
        else:
            distance_each_batch = precedence_constrained_routing(
                item_list,
                aisle,
                each_aisle_item,
                rack_x,
                rack_y,
                aisle_x,
                enter_aisle,
                distance_y,
            )
        list_distance_batch.append(distance_each_batch)

    number_item_each_batch = [len(batch) for batch in list_index_item_in_batch]

    setup_time = int(settings.get('setup_time', 0))
    picking_and_searching_time = int(settings.get('picking_and_searching_time', 0))
    speed_picker = int(settings.get('speed_picker', 1))
    sorting_time = int(settings.get('sorting_time', 0))

    process_time_batch: List[float] = []
    for distance, batch_size in zip(list_distance_batch, number_item_each_batch):
        process_time = calculate_completion_time(
            distance,
            batch_size,
            setup_time,
            picking_and_searching_time,
            speed_picker,
            sorting_time,
        )
        process_time_batch.append(float(format(process_time, '.3f')))

    for j, process_time in enumerate(process_time_batch, start=1):
        df_item_poor_new_PSO.loc[df_item_poor_new_PSO['batch'] == j, 'process_time'] = process_time

    return BatchPreparation(df_item_poor_new_PSO, list_batching_item_PSO, list_index_item_in_batch)


def evaluate_all_sols(
    new_pop_sol,
    df_item_poor_batch,
    heavy_item_set,
    name_path_input,
    scheduler: Optional[Any] = None,
):
    path_folder = os.path.join(name_path_input, f'setting_parameter_{name_path_input}.csv')
    df_setting_parameter = pd.read_csv(path_folder)

    settings: Dict[str, Any] = {}
    with open(path_folder) as file_obj:
        reader_obj = csv.reader(file_obj)
        for row in reader_obj:
            if not row:
                continue
            settings[row[0]] = _parse_setting_value(row[1]) if len(row) > 1 else None

    preparation = _prepare_batches(
        new_pop_sol,
        df_item_poor_batch,
        heavy_item_set,
        settings,
        name_path_input,
    )

    num_picker = int(settings.get('num_picker', 1))

    if scheduler is None:
        list_picker, list_tardiness_each_order, total_tardiness, df_item_poor_AS = ESDR_algorithms(
            preparation.dataframe,
            num_picker,
        )
    else:
        list_picker, list_tardiness_each_order, total_tardiness, df_item_poor_AS = scheduler(
            preparation.dataframe,
            num_picker,
        )

    return list_picker, total_tardiness, preparation.index_items_in_batch
