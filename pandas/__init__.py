"""A very small subset of pandas functionality used for the WaOA/AHA demos.
This module is **not** a drop-in replacement for pandas.  It only implements
just enough features for the accompanying optimisation scripts to run in the
execution environment where the real pandas package is unavailable.
"""

from __future__ import annotations

import csv
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union

Number = Union[int, float]


class _Index:
    def __init__(self, labels: Optional[Iterable[Any]] = None) -> None:
        self._labels: List[Any] = list(labels) if labels is not None else []

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._labels)

    def __iter__(self) -> Iterator[Any]:  # pragma: no cover - trivial
        return iter(self._labels)

    def __getitem__(self, item: Union[int, slice]) -> Any:
        return self._labels[item]

    def __setitem__(self, key: int, value: Any) -> None:
        self._labels[key] = value

    def __contains__(self, item: Any) -> bool:  # pragma: no cover - trivial
        return item in self._labels

    def append(self, value: Any) -> None:
        self._labels.append(value)

    def extend(self, values: Iterable[Any]) -> None:
        self._labels.extend(values)

    def remove(self, value: Any) -> None:
        self._labels.remove(value)

    def tolist(self) -> List[Any]:  # pragma: no cover - trivial
        return list(self._labels)

    @property
    def values(self) -> List[Any]:  # pragma: no cover - trivial
        return list(self._labels)

    def copy(self) -> "_Index":
        return _Index(self._labels)


class Series:
    def __init__(self, data: Iterable[Any], index: Optional[Iterable[Any]] = None, name: Optional[str] = None) -> None:
        self._data = list(data)
        if index is None:
            self.index = _Index(range(len(self._data)))
        else:
            self.index = _Index(index)
        self.name = name

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._data)

    def __iter__(self) -> Iterator[Any]:  # pragma: no cover - trivial
        return iter(self._data)

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        return f"Series({self._data})"

    def copy(self) -> "Series":
        return Series(self._data, self.index, self.name)

    def tolist(self) -> List[Any]:  # pragma: no cover - trivial
        return list(self._data)

    @property
    def values(self) -> "_SeriesValues":
        return _SeriesValues(self)

    def sum(self) -> Number:
        return sum(self._data)

    def max(self) -> Any:
        return max(self._data)

    def min(self) -> Any:
        return min(self._data)

    def astype(self, dtype: Any) -> "Series":
        if dtype is int:
            return Series([int(v) for v in self._data], self.index, self.name)
        if dtype is float:
            return Series([float(v) for v in self._data], self.index, self.name)
        raise TypeError(f"Unsupported dtype conversion: {dtype}")

    def isin(self, values: Iterable[Any]) -> "Series":
        value_set = set(values)
        return Series([val in value_set for val in self._data], self.index)

    def __getitem__(self, key: Union[int, slice, List[int], List[bool], "Series"]):
        if isinstance(key, Series):
            key = key.tolist()
        if isinstance(key, str):
            labels = self.index.tolist()
            if key not in labels:
                raise KeyError(key)
            position = labels.index(key)
            return self._data[position]
        if isinstance(key, int) and key in self.index.tolist():
            position = self.index.tolist().index(key)
            return self._data[position]
        if isinstance(key, list):
            if key and isinstance(key[0], bool):
                filtered = [val for val, flag in zip(self._data, key) if flag]
                filtered_index = [idx for idx, flag in zip(self.index, key) if flag]
                return Series(filtered, filtered_index, self.name)
            else:
                filtered = [self._data[i] for i in key]
                filtered_index = [self.index[i] for i in key]
                return Series(filtered, filtered_index, self.name)
        if isinstance(key, slice):
            return Series(self._data[key], self.index[key], self.name)
        return self._data[key]

    def __setitem__(self, key: int, value: Any) -> None:
        self._data[key] = value

    def _binary_op(self, other: Any, op) -> "Series":
        if isinstance(other, Series):
            other_values = other._data
        else:
            other_values = [other] * len(self._data)
        return Series([op(a, b) for a, b in zip(self._data, other_values)], self.index, self.name)

    def _compare(self, other: Any, comparator) -> "Series":
        if isinstance(other, Series):
            other_values = other._data
        else:
            other_values = [other] * len(self._data)
        results = []
        for a, b in zip(self._data, other_values):
            if a is None or b is None:
                results.append(False)
            else:
                results.append(comparator(a, b))
        return Series(results, self.index, self.name)

    def __add__(self, other: Any) -> "Series":
        return self._binary_op(other, lambda a, b: a + b)

    def __radd__(self, other: Any) -> "Series":
        return self.__add__(other)

    def __truediv__(self, other: Number) -> "Series":
        return Series([val / other for val in self._data], self.index, self.name)

    def __eq__(self, other: Any) -> "Series":
        return self._binary_op(other, lambda a, b: a == b)

    def __ne__(self, other: Any) -> "Series":
        return self._binary_op(other, lambda a, b: a != b)

    def __lt__(self, other: Any) -> "Series":
        return self._compare(other, lambda a, b: a < b)

    def __le__(self, other: Any) -> "Series":
        return self._compare(other, lambda a, b: a <= b)

    def __gt__(self, other: Any) -> "Series":
        return self._compare(other, lambda a, b: a > b)

    def __ge__(self, other: Any) -> "Series":
        return self._compare(other, lambda a, b: a >= b)


class _SeriesValues:
    def __init__(self, series: Series) -> None:
        self._series = series

    def __iter__(self) -> Iterator[Any]:  # pragma: no cover - trivial
        return iter(self._series._data)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._series._data)

    def __getitem__(self, item: Union[int, slice]) -> Any:
        if isinstance(item, slice):
            subset = Series(self._series._data[item], self._series.index[item], self._series.name)
            return _SeriesValues(subset)
        return self._series._data[item]

    def __float__(self) -> float:
        if len(self._series._data) != 1:
            raise TypeError("Only length-1 values can be converted to float")
        return float(self._series._data[0])

    def __int__(self) -> int:
        if len(self._series._data) != 1:
            raise TypeError("Only length-1 values can be converted to int")
        return int(self._series._data[0])

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        return repr(self._series._data)

    def tolist(self) -> List[Any]:  # pragma: no cover - trivial
        return list(self._series._data)

    def sum(self) -> Number:
        return sum(self._series._data)

    def max(self) -> Any:
        return max(self._series._data)

    def min(self) -> Any:
        return min(self._series._data)

    def __add__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a + b)

    def __radd__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: b + a)

    def __sub__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a - b)

    def __rsub__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: b - a)

    def __mul__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a * b)

    def __truediv__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a / b)

    def __rtruediv__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: b / a)

    def __lt__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a < b)

    def __le__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a <= b)

    def __gt__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a > b)

    def __ge__(self, other: Any) -> Any:
        return self._binary_op(other, lambda a, b: a >= b)

    def _binary_op(self, other: Any, op) -> Any:
        if isinstance(other, _SeriesValues):
            other_iter = other._series._data
        elif isinstance(other, Series):
            other_iter = other._data
        elif isinstance(other, list):
            other_iter = other
        else:
            other_iter = [other] * len(self._series._data)
        result = [op(a, b) for a, b in zip(self._series._data, other_iter)]
        if len(result) == 1:
            return result[0]
        return result


class DataFrame:
    def __init__(self, data: Optional[Iterable[Any]] = None, columns: Optional[Sequence[Any]] = None,
                 index: Optional[Iterable[Any]] = None) -> None:
        if data is None:
            self._data: List[Dict[Any, Any]] = []
            self.columns = _Index(columns or [])
            self.index = _Index(index or [])
            return

        if isinstance(data, DataFrame):
            self._data = [row.copy() for row in data._data]
            self.columns = data.columns.copy()
            self.index = data.index.copy()
            return

        rows: List[Dict[Any, Any]] = []
        if isinstance(data, list) and data and isinstance(data[0], dict):
            rows = [row.copy() for row in data]
            if columns is not None:
                col_names = list(columns)
            else:
                col_names = []
                for row in rows:
                    for key in row:
                        if key not in col_names:
                            col_names.append(key)
        elif isinstance(data, list) and data and not isinstance(data[0], dict):
            col_names = list(columns) if columns is not None else list(range(len(data[0])))
            for values in data:
                row = {col_names[i]: values[i] for i in range(len(col_names))}
                rows.append(row)
        else:
            rows = []
            col_names = list(columns or [])

        self._data = rows
        self.columns = _Index(col_names)
        if index is None:
            self.index = _Index(range(len(rows)))
        else:
            self.index = _Index(index)

    def __len__(self) -> int:  # pragma: no cover - trivial
        return len(self._data)

    def __repr__(self) -> str:  # pragma: no cover - debugging helper
        return f"DataFrame({self._data})"

    @property
    def shape(self) -> Tuple[int, int]:  # pragma: no cover - trivial
        return len(self._data), len(self.columns)

    def copy(self) -> "DataFrame":
        return DataFrame(self)

    def to_dicts(self) -> List[Dict[Any, Any]]:  # pragma: no cover - debugging helper
        return [row.copy() for row in self._data]

    def _ensure_column(self, column: Any) -> None:
        if column in self.columns:
            return
        self.columns.append(column)
        for row in self._data:
            row[column] = None

    def _column_values(self, column: Any) -> List[Any]:
        self._ensure_column(column)
        return [row.get(column) for row in self._data]

    def __getitem__(self, key: Union[str, int, List[Any], Series]) -> Union[Series, "DataFrame"]:
        if isinstance(key, Series) and all(isinstance(val, bool) for val in key._data):
            rows = [row.copy() for row, flag in zip(self._data, key._data) if flag]
            index = [idx for idx, flag in zip(self.index, key._data) if flag]
            return DataFrame(rows, columns=self.columns, index=index)
        if isinstance(key, Series):
            key = key.tolist()
        if isinstance(key, list):
            if key and isinstance(key[0], bool):
                rows = [row.copy() for row, flag in zip(self._data, key) if flag]
                index = [idx for idx, flag in zip(self.index, key) if flag]
                return DataFrame(rows, columns=self.columns, index=index)
            return DataFrame([{col: row.get(col) for col in key} for row in self._data], columns=key, index=self.index)
        values = self._column_values(key)
        return Series(values, self.index, name=key)

    def __setitem__(self, key: Any, value: Any) -> None:
        self._ensure_column(key)
        if isinstance(value, Series):
            value_list = value.tolist()
        elif isinstance(value, list):
            value_list = value
        else:
            value_list = [value] * len(self._data)
        for row, val in zip(self._data, value_list):
            row[key] = val

    @property
    def loc(self) -> "_LocIndexer":
        return _LocIndexer(self)

    @property
    def iloc(self) -> "_ILocIndexer":
        return _ILocIndexer(self)

    def drop(self, labels: Iterable[Any], axis: int = 0) -> "DataFrame":
        if axis == 0:
            labels_set = set(labels)
            new_rows = []
            new_index = []
            for idx_label, row in zip(self.index, self._data):
                if idx_label not in labels_set:
                    new_rows.append(row.copy())
                    new_index.append(idx_label)
            result = DataFrame(new_rows, columns=self.columns, index=new_index)
            return result
        elif axis == 1:
            drop_set = set(labels)
            new_columns = [col for col in self.columns if col not in drop_set]
            new_rows = []
            for row in self._data:
                new_rows.append({col: row.get(col) for col in new_columns})
            return DataFrame(new_rows, columns=new_columns, index=self.index)
        else:
            raise ValueError("axis must be 0 or 1")

    def reset_index(self, inplace: bool = False) -> Optional["DataFrame"]:
        old_index = self.index.tolist()
        new_rows = []
        for new_idx, (old_idx, row) in enumerate(zip(old_index, self._data)):
            new_row = row.copy()
            new_row['index'] = old_idx
            new_rows.append(new_row)
        result = DataFrame(new_rows, columns=list(self.columns) + ['index'], index=range(len(new_rows)))
        if inplace:
            self._data = result._data
            self.columns = result.columns
            self.index = result.index
            return None
        return result

    def sort_values(self, by: Union[str, List[str]], ascending: Union[bool, List[bool]] = True, inplace: bool = False) -> "DataFrame":
        if isinstance(by, list):
            sort_keys = by
        else:
            sort_keys = [by]
        if isinstance(ascending, list):
            ascendings = ascending
        else:
            ascendings = [ascending] * len(sort_keys)

        new_rows = [row.copy() for row in self._data]
        for col, asc in reversed(list(zip(sort_keys, ascendings))):
            new_rows.sort(key=lambda row: row.get(col), reverse=not asc)
        result = DataFrame(new_rows, columns=self.columns, index=self.index)
        if inplace:
            self._data = result._data
            return self
        return result

    def assign(self, **kwargs: Any) -> "DataFrame":
        new_rows = [row.copy() for row in self._data]
        new_columns = list(self.columns)
        for column, value in kwargs.items():
            if isinstance(value, list):
                values = value
            else:
                values = [value] * len(new_rows)
            for row, val in zip(new_rows, values):
                row[column] = val
            if column not in new_columns:
                new_columns.append(column)
        return DataFrame(new_rows, columns=new_columns, index=self.index)

    def astype(self, dtype: Any) -> "DataFrame":
        if isinstance(dtype, dict):
            new_rows = []
            for row in self._data:
                new_row = row.copy()
                for column, target in dtype.items():
                    if column in new_row and new_row[column] is not None:
                        new_row[column] = target(new_row[column])
                new_rows.append(new_row)
            return DataFrame(new_rows, columns=self.columns, index=self.index)

        if dtype is int:
            new_rows = []
            for row in self._data:
                new_rows.append({key: int(value) if value is not None else value for key, value in row.items()})
            return DataFrame(new_rows, columns=self.columns, index=self.index)

        if dtype is float:
            new_rows = []
            for row in self._data:
                new_rows.append({key: float(value) if value is not None else value for key, value in row.items()})
            return DataFrame(new_rows, columns=self.columns, index=self.index)

        raise TypeError("Unsupported dtype conversion")

    def __getattr__(self, name: str) -> Any:
        if name == 'values':
            return [row.copy() for row in self._data]
        raise AttributeError(name)


class _LocIndexer:
    def __init__(self, df: DataFrame) -> None:
        self.df = df

    def _resolve_rows(self, row_sel: Any) -> Tuple[List[int], List[Any]]:
        if isinstance(row_sel, Series):
            mask = row_sel.tolist()
            indices = [i for i, flag in enumerate(mask) if flag]
            labels = [self.df.index[i] for i in indices]
            return indices, labels
        if isinstance(row_sel, list):
            if row_sel and isinstance(row_sel[0], bool):
                indices = [i for i, flag in enumerate(row_sel) if flag]
                labels = [self.df.index[i] for i in indices]
                return indices, labels
            labels = row_sel
            indices = [self.df.index.tolist().index(label) for label in labels]
            return indices, labels
        if isinstance(row_sel, slice):
            start = row_sel.start or 0
            stop = row_sel.stop if row_sel.stop is not None else len(self.df.index)
            indices = list(range(start, stop))
            labels = [self.df.index[i] for i in indices]
            return indices, labels
        if row_sel is None:
            indices = list(range(len(self.df.index)))
            labels = self.df.index.tolist()
            return indices, labels
        # single label or append
        if row_sel in self.df.index.tolist():
            position = self.df.index.tolist().index(row_sel)
            return [position], [row_sel]
        if isinstance(row_sel, int) and row_sel == len(self.df.index):
            # append a new row with this label
            self.df.index.append(row_sel)
            self.df._data.append({col: None for col in self.df.columns})
            return [row_sel], [row_sel]
        # assume numeric positional index
        if isinstance(row_sel, int) and 0 <= row_sel < len(self.df.index):
            return [row_sel], [self.df.index[row_sel]]
        raise KeyError(f"Row selector {row_sel} not found")

    def __getitem__(self, key: Any) -> Union[Series, DataFrame]:
        if isinstance(key, tuple):
            row_sel, col_sel = key
        else:
            row_sel, col_sel = key, None
        row_indices, row_labels = self._resolve_rows(row_sel)

        if col_sel is None:
            rows = [self.df._data[i].copy() for i in row_indices]
            cols = list(self.df.columns)
            return DataFrame(rows, columns=cols, index=row_labels)

        if isinstance(col_sel, list):
            cols = col_sel
        elif isinstance(col_sel, slice):
            cols = self.df.columns[col_sel]
        else:
            cols = [col_sel]

        data = []
        for i in row_indices:
            row = self.df._data[i]
            data.append({col: row.get(col) for col in cols})
        result_df = DataFrame(data, columns=cols, index=row_labels)
        if len(cols) == 1:
            return result_df[cols[0]]
        return result_df

    def __setitem__(self, key: Any, value: Any) -> None:
        if isinstance(key, tuple):
            row_sel, col_sel = key
        else:
            row_sel, col_sel = key, None

        row_indices, row_labels = self._resolve_rows(row_sel)

        if col_sel is None:
            if isinstance(value, Series):
                row_values_list = [value.tolist()]
            elif isinstance(value, list) and row_indices and isinstance(value[0], list):
                row_values_list = value
            elif isinstance(value, list):
                row_values_list = [value]
            else:
                raise ValueError("Expected list or Series for row assignment")

            if len(row_values_list) == 1 and len(row_indices) > 1:
                row_values_list = row_values_list * len(row_indices)

            for row_idx, values in zip(row_indices, row_values_list):
                if len(values) != len(self.df.columns):
                    raise ValueError("Row length mismatch")
                for col, val in zip(self.df.columns, values):
                    self.df._data[row_idx][col] = val
            return

        if isinstance(col_sel, list):
            columns = col_sel
        elif isinstance(col_sel, slice):
            columns = self.df.columns[col_sel]
        else:
            columns = [col_sel]

        for column in columns:
            self.df._ensure_column(column)

        if isinstance(value, Series):
            values_iter = value.tolist()
            if len(columns) != 1:
                raise ValueError("Series assignment only supported for single column")
            for row_idx, val in zip(row_indices, values_iter):
                self.df._data[row_idx][columns[0]] = val
            return

        if isinstance(value, list) and len(value) == len(row_indices):
            value_rows = value
        else:
            value_rows = [value] * len(row_indices)

        for row_idx, row_value in zip(row_indices, value_rows):
            if len(columns) == 1:
                self.df._data[row_idx][columns[0]] = row_value
            else:
                for col, val in zip(columns, row_value):
                    self.df._data[row_idx][col] = val


class _ILocIndexer:
    def __init__(self, df: DataFrame) -> None:
        self.df = df

    def __getitem__(self, key: Any) -> Union[Series, DataFrame]:
        if isinstance(key, tuple):
            row_sel, col_sel = key
        else:
            row_sel, col_sel = key, None

        is_single_label = False
        if isinstance(row_sel, list):
            rows = [self.df._data[i].copy() for i in row_sel]
            index = [self.df.index[i] for i in row_sel]
        elif isinstance(row_sel, slice):
            indices = list(range(*row_sel.indices(len(self.df._data))))
            rows = [self.df._data[i].copy() for i in indices]
            index = [self.df.index[i] for i in indices]
        elif isinstance(row_sel, int):
            rows = [self.df._data[row_sel].copy()]
            index = [self.df.index[row_sel]]
            is_single_label = True
        else:
            raise KeyError(f"Unsupported iloc selector: {row_sel}")

        if col_sel is None:
            if len(rows) == 1 and is_single_label:
                series_data = [rows[0].get(col) for col in self.df.columns]
                return Series(series_data, index=self.df.columns, name=index[0])
            return DataFrame(rows, columns=self.df.columns, index=index)

        if isinstance(col_sel, list):
            cols = col_sel
        elif isinstance(col_sel, slice):
            cols = self.df.columns[col_sel]
        else:
            cols = [col_sel]

        data = [{col: row.get(col) for col in cols} for row in rows]
        result = DataFrame(data, columns=cols, index=index)
        if len(cols) == 1:
            return result[cols[0]]
        return result


def read_csv(path: str, header: Optional[int] = 0) -> DataFrame:
    with open(path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)
    if header is None:
        columns = list(range(len(rows[0]))) if rows else []
        data_rows = rows
    else:
        columns = rows[0] if rows else []
        data_rows = rows[1:]
    dict_rows = [
        {columns[i]: _convert_value(row[i]) for i in range(len(columns)) if i < len(row)}
        for row in data_rows
    ]
    return DataFrame(dict_rows, columns=columns)


def _convert_value(value: str) -> Any:
    if value == '':
        return 0
    try:
        if '.' in value:
            return float(value)
        return int(value)
    except ValueError:
        return value


def concat(dfs: List[DataFrame], ignore_index: bool = False) -> DataFrame:
    if not dfs:
        return DataFrame()
    columns = list(dfs[0].columns)
    rows: List[Dict[Any, Any]] = []
    index: List[Any] = []
    for df in dfs:
        for row in df._data:
            rows.append(row.copy())
        index.extend(df.index)
    if ignore_index:
        index = list(range(len(rows)))
    return DataFrame(rows, columns=columns, index=index)


def DataFrame_from_records(records: List[Dict[Any, Any]]) -> DataFrame:  # pragma: no cover - helper
    return DataFrame(records)


def set_option(*args: Any, **kwargs: Any) -> None:  # pragma: no cover - no-op for compatibility
    return None


__all__ = [
    "DataFrame",
    "Series",
    "concat",
    "read_csv",
    "set_option",
]
