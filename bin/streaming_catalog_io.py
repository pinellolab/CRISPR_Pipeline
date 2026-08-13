#!/usr/bin/env python3
"""Streaming fast paths for enriched Parquet catalog tables."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


def _is_parquet(path: str | Path) -> bool:
    return str(path).lower().endswith((".parquet", ".pq"))


def _normalize_join_keys(local_scan, global_scan, join_columns, pl):
    """Cast equivalent dictionary/string and numeric join keys identically.

    Parquet written from pandas categoricals is scanned by Polars as
    ``Categorical``, whereas the same identifiers in Polars-written output are
    ``String``. Polars deliberately refuses to guess this cast during a join.
    Normalize only the key columns and preserve every metric's native dtype.
    """

    local_schema = local_scan.collect_schema()
    global_schema = global_scan.collect_schema()
    string_types = (pl.String, pl.Categorical, pl.Enum)
    local_exprs = []
    global_exprs = []
    for column in join_columns:
        local_dtype = local_schema[column]
        global_dtype = global_schema[column]
        if local_dtype in string_types or global_dtype in string_types:
            common_dtype = pl.String
        elif local_dtype.is_integer() and global_dtype.is_integer():
            common_dtype = pl.Int64
        elif (
            (local_dtype.is_integer() or local_dtype.is_float())
            and (global_dtype.is_integer() or global_dtype.is_float())
        ):
            common_dtype = pl.Float64
        elif local_dtype == global_dtype:
            continue
        else:
            raise TypeError(
                f"Cannot normalize catalog join column {column!r}: "
                f"local={local_dtype}, global={global_dtype}"
            )
        local_exprs.append(pl.col(column).cast(common_dtype))
        global_exprs.append(pl.col(column).cast(common_dtype))

    if local_exprs:
        local_scan = local_scan.with_columns(local_exprs)
    if global_exprs:
        global_scan = global_scan.with_columns(global_exprs)
    return local_scan, global_scan


def try_write_enriched_parquet_catalog(
    *,
    local_path: str | Path,
    global_path: str | Path,
    output_path: str | Path,
    join_columns: Sequence[str],
    local_metric_columns: Sequence[str],
    global_required_columns: Sequence[str],
    output_columns: Sequence[str],
    sort_columns: Sequence[str],
) -> bool:
    """Write an enriched catalog without materializing all pairs in pandas."""

    if not all(_is_parquet(path) for path in (local_path, global_path, output_path)):
        return False

    try:
        import polars as pl
    except ImportError:
        return False

    local_scan = pl.scan_parquet(local_path)
    global_scan = pl.scan_parquet(global_path)
    local_schema = set(local_scan.collect_schema().names())
    global_schema = set(global_scan.collect_schema().names())
    local_required = set(join_columns) | set(local_metric_columns)
    global_required = set(join_columns) | set(global_required_columns)
    if not local_required.issubset(local_schema):
        return False
    if not global_required.issubset(global_schema):
        return False

    local_scan, global_scan = _normalize_join_keys(
        local_scan, global_scan, join_columns, pl
    )

    local_keys = local_scan.select(join_columns)
    duplicate_key = (
        local_keys.group_by(join_columns)
        .len()
        .filter(pl.col("len") > 1)
        .limit(1)
        .collect(engine="streaming")
    )
    if duplicate_key.height:
        return False

    local_only = (
        local_keys.join(
            global_scan.select(join_columns),
            on=join_columns,
            how="anti",
        )
        .limit(1)
        .collect(engine="streaming")
    )
    if local_only.height:
        return False

    local_metrics = local_scan.select(list(join_columns) + list(local_metric_columns))
    catalog = global_scan.join(local_metrics, on=join_columns, how="left")
    catalog = catalog.select(output_columns).sort(
        sort_columns,
        nulls_last=True,
        maintain_order=True,
    )
    catalog.sink_parquet(
        output_path,
        compression="zstd",
        maintain_order=True,
        engine="streaming",
    )
    print(
        "Wrote catalog with the streaming enriched-Parquet fast path: "
        f"{output_path}",
        flush=True,
    )
    return True
