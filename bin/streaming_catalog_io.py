#!/usr/bin/env python3
"""Streaming fast paths for enriched Parquet catalog tables."""

from __future__ import annotations

from pathlib import Path
from typing import Sequence


def _is_parquet(path: str | Path) -> bool:
    return str(path).lower().endswith((".parquet", ".pq"))


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
