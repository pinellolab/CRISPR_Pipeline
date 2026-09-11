#!/usr/bin/env python
"""Patch a MuData file's top-level .uns without rewriting its assay matrices."""

from __future__ import annotations

from pathlib import Path
import shutil
from typing import Iterable

import h5py
import numpy as np

try:
    from anndata.io import write_elem
except ImportError:  # older anndata
    from anndata.experimental import write_elem


def write_uns_patch(
    input_path: str | Path,
    output_path: str | Path,
    updates: dict | None = None,
    deletes: Iterable[str] = (),
    aliases: dict[str, str] | None = None,
) -> None:
    """Copy ``input_path`` to ``output_path``, then add/overwrite ``updates``
    and remove ``deletes`` from the file's top-level /uns group.

    Never touches /mod/*/X or /mod/*/layers, so this is far cheaper than a
    full MuData read + write when the only change is to .uns -- the copy is
    a raw byte copy (no re-parsing/re-chunking/re-compression), and the /uns
    patch only re-serializes whatever was already in .uns plus `updates`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if input_path.resolve() != output_path.resolve():
        shutil.copy(input_path, output_path)

    with h5py.File(output_path, "r+") as f:
        # Patch individual children rather than reading and rewriting the whole
        # /uns mapping.  In global screens a single result table can contain
        # >100M rows; reconstructing the full mapping multiplies memory use and
        # makes each subsequent update rewrite every table already attached.
        if "uns" not in f:
            write_elem(f, "uns", {})
        uns = f["uns"]
        for key in deletes:
            if key in uns:
                del uns[key]
        for key, value in (updates or {}).items():
            if key in uns:
                del uns[key]
            write_elem(uns, key, value)
        # An alias is an HDF5 hard link, not a second copy: the two names point at
        # one set of datasets, so a consumer reading either sees the same table and
        # the file does not carry it twice. A screen-scale result table is gigabytes,
        # and the generic keys are exactly the same frames as the analysis-qualified
        # ones, so writing both cost several gigabytes of pure duplication.
        for alias, target in (aliases or {}).items():
            if target not in uns:
                raise KeyError(f"cannot alias {alias!r} to {target!r}: {target!r} is not in /uns")
            if alias in uns:
                del uns[alias]
            uns[alias] = uns[target]


def _categorical_code_dtype(category_count: int) -> np.dtype:
    if category_count <= np.iinfo(np.int8).max + 1:
        return np.dtype("int8")
    if category_count <= np.iinfo(np.int16).max + 1:
        return np.dtype("int16")
    if category_count <= np.iinfo(np.int32).max + 1:
        return np.dtype("int32")
    return np.dtype("int64")


def write_parquet_dataframe_to_uns(
    h5mu_path: str | Path,
    key: str,
    parquet_path: str | Path,
) -> None:
    """Stream a Parquet table directly into an AnnData dataframe in ``/uns``.

    Columns are decoded independently with Polars. The complete table never
    exists as pandas or as a single in-memory Arrow table; peak memory is
    bounded by the largest single column. A temporary HDF5 group is atomically
    renamed on success, so an interrupted write cannot replace an existing
    result table with a partial one.
    """

    import gc
    import polars as pl

    h5mu_path = Path(h5mu_path)
    parquet_path = Path(parquet_path)
    scan = pl.scan_parquet(parquet_path)
    schema = scan.collect_schema()
    column_names = schema.names()
    row_count = (
        scan.select(pl.len().alias("rows"))
        .collect(engine="streaming")
        .item()
    )
    temporary_key = f"__{key}_streaming_tmp"
    string_dtype = h5py.string_dtype(encoding="utf-8")
    polars_signed_types = {
        1: pl.Int8,
        2: pl.Int16,
        4: pl.Int32,
        8: pl.Int64,
    }

    with h5py.File(h5mu_path, "r+") as handle:
        if "uns" not in handle:
            write_elem(handle, "uns", {})
        uns = handle["uns"]
        if temporary_key in uns:
            del uns[temporary_key]
        frame = uns.create_group(temporary_key)
        frame.attrs["encoding-type"] = "dataframe"
        frame.attrs["encoding-version"] = "0.2.0"
        frame.attrs["_index"] = "_index"
        frame.attrs.create(
            "column-order",
            np.asarray(column_names, dtype=object),
            dtype=string_dtype,
        )

        index = frame.create_dataset(
            "_index", shape=(row_count,), dtype=np.int64, chunks=True
        )
        index.attrs["encoding-type"] = "array"
        index.attrs["encoding-version"] = "0.2.0"
        index_chunk_rows = 1_000_000
        for start in range(0, row_count, index_chunk_rows):
            stop = min(start + index_chunk_rows, row_count)
            index[start:stop] = np.arange(start, stop, dtype=np.int64)

        for column_index, name in enumerate(column_names, start=1):
            source_dtype = schema[name]
            if source_dtype in (pl.String, pl.Categorical, pl.Enum):
                series = (
                    scan.select(pl.col(name).cast(pl.Categorical))
                    .collect(engine="streaming")[name]
                )
                categories_values = series.cat.get_categories().to_list()
                code_dtype = _categorical_code_dtype(len(categories_values))
                polars_code_dtype = polars_signed_types[code_dtype.itemsize]
                codes_values = (
                    series.to_physical()
                    .cast(polars_code_dtype)
                    .fill_null(-1)
                    .to_numpy()
                )
                group = frame.create_group(name)
                group.attrs["encoding-type"] = "categorical"
                group.attrs["encoding-version"] = "0.2.0"
                group.attrs["ordered"] = False
                codes = group.create_dataset(
                    "codes", data=codes_values, chunks=True
                )
                codes.attrs["encoding-type"] = "array"
                codes.attrs["encoding-version"] = "0.2.0"
                categories = group.create_dataset(
                    "categories",
                    data=np.asarray(categories_values, dtype=object),
                    dtype=string_dtype,
                )
                categories.attrs["encoding-type"] = "string-array"
                categories.attrs["encoding-version"] = "0.2.0"
                del series, categories_values, codes_values
            else:
                series = scan.select(name).collect(engine="streaming")[name]
                has_nulls = series.null_count() > 0
                if source_dtype.is_integer() or source_dtype == pl.Boolean:
                    if has_nulls:
                        group = frame.create_group(name)
                        encoding_type = (
                            "nullable-boolean"
                            if source_dtype == pl.Boolean
                            else "nullable-integer"
                        )
                        group.attrs["encoding-type"] = encoding_type
                        group.attrs["encoding-version"] = "0.1.0"
                        values_array = series.fill_null(
                            False if source_dtype == pl.Boolean else 0
                        ).to_numpy()
                        mask_array = series.is_null().to_numpy()
                        values = group.create_dataset(
                            "values", data=values_array, chunks=True
                        )
                        mask = group.create_dataset(
                            "mask", data=mask_array, chunks=True
                        )
                        for dataset in (values, mask):
                            dataset.attrs["encoding-type"] = "array"
                            dataset.attrs["encoding-version"] = "0.2.0"
                        del values_array, mask_array
                    else:
                        values_array = series.to_numpy()
                        dataset = frame.create_dataset(
                            name, data=values_array, chunks=True
                        )
                        dataset.attrs["encoding-type"] = "array"
                        dataset.attrs["encoding-version"] = "0.2.0"
                        del values_array
                elif source_dtype.is_float():
                    values_array = series.fill_null(float("nan")).to_numpy()
                    dataset = frame.create_dataset(
                        name, data=values_array, chunks=True
                    )
                    dataset.attrs["encoding-type"] = "array"
                    dataset.attrs["encoding-version"] = "0.2.0"
                    del values_array
                else:
                    raise TypeError(
                        f"Unsupported Polars type for H5MU dataframe column "
                        f"{name!r}: {source_dtype}"
                    )
                del series

            gc.collect()
            handle.flush()
            print(
                f"Embedded {key}: column {column_index}/{len(column_names)} "
                f"({name}; {row_count} rows)",
                flush=True,
            )

        if key in uns:
            del uns[key]
        uns.move(temporary_key, key)
