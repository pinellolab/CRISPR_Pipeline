#!/usr/bin/env python3
"""Streaming fast paths for enriched Parquet catalog tables."""

from __future__ import annotations

from polars_compat import lazy_schema, JOIN_ORDER_LEFT, UNIQUE_ORDER, SORT_ORDER, SINK_ORDER, SINK_ENGINE, COLLECT_ENGINE

from pathlib import Path
from typing import Mapping, Sequence


def _is_parquet(path: str | Path) -> bool:
    return str(path).lower().endswith((".parquet", ".pq"))


def _normalize_join_keys(local_scan, global_scan, join_columns, pl):
    """Cast equivalent dictionary/string and numeric join keys identically.

    Parquet written from pandas categoricals is scanned by Polars as
    ``Categorical``, whereas the same identifiers in Polars-written output are
    ``String``. Polars deliberately refuses to guess this cast during a join.
    Normalize only the key columns and preserve every metric's native dtype.
    """

    local_schema = lazy_schema(local_scan)
    global_schema = lazy_schema(global_scan)
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


def try_write_enriched_parquet_catalog(**kwargs) -> bool:
    """Try the Polars fast path, and say so plainly when it cannot be taken.

    Every ``return False`` here means "the caller should build this catalog the
    ordinary way", and the pandas path is the reference implementation, so
    falling back is always correct and only ever slower. A Polars limitation on
    the input -- the Gasperini run met ``not implemented: reading dictionaries
    of type (Int32, Null)``, an all-null column -- is that same situation
    arriving as an exception rather than as a schema check, so it is handled the
    same way instead of failing the process.
    """
    try:
        return _write_enriched_parquet_catalog(**kwargs)
    except Exception as exc:  # noqa: BLE001 - any Polars limitation falls back
        print(
            "Streaming enriched-Parquet catalog unavailable, using the pandas "
            f"path instead: {type(exc).__name__}: {exc}",
            flush=True,
        )
        return False


def _write_enriched_parquet_catalog(
    *,
    local_path: str | Path,
    global_path: str | Path,
    output_path: str | Path,
    join_columns: Sequence[str],
    local_metric_columns: Sequence[str],
    global_required_columns: Sequence[str],
    output_columns: Sequence[str],
    sort_columns: Sequence[str],
    local_alias_columns: Mapping[str, str] | None = None,
    derived_neglog10: Mapping[str, str] | None = None,
    require_fillable_q: Sequence[tuple[str, str]] = (),
    pvalue_floor: float = 1e-300,
) -> bool:
    """Write an enriched catalog without materializing all pairs in pandas.

    ``local_alias_columns`` maps a column in the local table to the name it takes
    in the output, which is how the ``perturbo_cis_*`` block is produced: the
    pandas path renames the local table's ``perturbo_*`` columns, so those names
    exist in neither input file and the fast path has to perform the same rename
    rather than select a column that was never there. A mapped column the local
    table does not have is filled with nulls, matching the pandas path's
    ``np.nan`` fill for a run with no local PerTurbo results.

    ``derived_neglog10`` maps an output column to the p-value column it is
    computed from, again mirroring pandas rather than carrying a precomputed
    ``negLog10p`` across, so the floor applied here is the caller's.
    """

    if not all(_is_parquet(path) for path in (local_path, global_path, output_path)):
        return False

    try:
        import polars as pl
    except ImportError:
        return False

    local_scan = pl.scan_parquet(local_path)
    global_scan = pl.scan_parquet(global_path)
    local_schema = set(lazy_schema(local_scan).names())
    global_schema = set(lazy_schema(global_scan).names())
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
        .collect(**COLLECT_ENGINE)
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
        .collect(**COLLECT_ENGINE)
    )
    if local_only.height:
        return False

    alias = dict(local_alias_columns or {})
    derived = dict(derived_neglog10 or {})
    present_alias = {src: dst for src, dst in alias.items() if src in local_schema}
    # A mapped column the local table lacks is filled with nulls -- except where a
    # derivation below can rebuild it from a column that is present, which must be
    # left alone for that derivation to see it as missing.
    missing_alias = [
        dst
        for src, dst in alias.items()
        if src not in local_schema and dst not in derived
    ]

    local_metrics = local_scan.select(
        list(join_columns) + list(local_metric_columns) + list(present_alias)
    )
    if present_alias:
        local_metrics = local_metrics.rename(present_alias)
    catalog = global_scan.join(local_metrics, on=join_columns, how="left")

    if missing_alias:
        catalog = catalog.with_columns(
            [pl.lit(None, dtype=pl.Float64).alias(name) for name in missing_alias]
        )

    # Derive only what nothing else produced. mergedResults computes the local
    # table's negLog10p with the very numpy helper the pandas catalog path uses,
    # on the same p-value and the same floor, so carrying that column across is
    # bit-identical to recomputing it, whereas Polars' log10 differs from numpy's
    # in the last bit. Derivation is the fallback for a table that lacks it.
    produced = set(lazy_schema(catalog).names())
    for out_column, p_column in derived.items():
        if out_column in produced:
            continue
        catalog = catalog.with_columns(
            pl.col(p_column)
            .cast(pl.Float64)
            .clip(lower_bound=pvalue_floor)
            .log10()
            .neg()
            .alias(out_column)
        )

    # The pandas path fills a missing q-value by running Benjamini-Hochberg over
    # the local family. mergedResults already did that before writing this table,
    # so there should be nothing left to fill; Benjamini-Hochberg needs a global
    # sort that does not belong in a streaming plan, so if any row still wants a
    # q-value, hand the table back to pandas rather than emit a null it would
    # have filled.
    for q_column, p_column in require_fillable_q:
        unfillable = (
            catalog.select(
                (pl.col(q_column).is_null() & pl.col(p_column).is_not_null()).any()
            )
            .collect(**COLLECT_ENGINE)
            .item()
        )
        if unfillable:
            return False

    catalog = catalog.select(output_columns).sort(
        sort_columns,
        nulls_last=True,
        **SORT_ORDER,
    )
    catalog.sink_parquet(
        output_path,
        compression="zstd",
        **SINK_ORDER,
        **SINK_ENGINE,
    )
    print(
        "Wrote catalog with the streaming enriched-Parquet fast path: "
        f"{output_path}",
        flush=True,
    )
    return True
