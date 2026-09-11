#!/usr/bin/env python3
"""Vectorized, streaming enrichment for large PerTurbo Parquet tables."""

from __future__ import annotations

from polars_compat import lazy_schema

from pathlib import Path
from typing import Literal

import mudata as mu
import pandas as pd

from analysis_output_formatting import (
    ELEMENT_ANNOTATION_COLUMNS,
    ELEMENT_COLUMNS,
    GUIDE_ANNOTATION_COLUMNS,
    PERTURBO_RESULT_COLUMNS,
    P_VALUE_FLOOR,
    _build_element_metadata,
    _build_gene_name_map,
    _build_guide_metadata,
    _element_key,
    _normalize_element_columns,
)


GUIDE_OUTPUT_COLUMNS = [
    "gene_id",
    "guide_id",
    *PERTURBO_RESULT_COLUMNS,
    *GUIDE_ANNOTATION_COLUMNS,
]
ELEMENT_OUTPUT_COLUMNS = [
    "gene_id",
    *ELEMENT_COLUMNS,
    *PERTURBO_RESULT_COLUMNS,
    *ELEMENT_ANNOTATION_COLUMNS,
]


def _normalize_metadata_text(frame: pd.DataFrame) -> pd.DataFrame:
    """Give Polars homogeneous strings without altering numeric columns."""

    result = frame.copy()
    for column in result.columns:
        series = result[column]
        if (
            pd.api.types.is_object_dtype(series.dtype)
            or pd.api.types.is_string_dtype(series.dtype)
            or isinstance(series.dtype, pd.CategoricalDtype)
        ):
            result[column] = series.astype("string").fillna("").astype(str)
    return result


def _guide_metadata(mdata) -> pd.DataFrame:
    return _normalize_metadata_text(_build_guide_metadata(mdata))


def _element_metadata(mdata) -> pd.DataFrame:
    metadata = _build_element_metadata(mdata)
    guide_var = _normalize_element_columns(mdata["guide"].var.copy())
    guide_var["_element_key"] = _element_key(guide_var)
    native_keys = guide_var[["_element_key", *ELEMENT_COLUMNS]].drop_duplicates(
        "_element_key", keep="first"
    )
    metadata = metadata.merge(native_keys, on="_element_key", how="left")
    return _normalize_metadata_text(metadata.drop(columns="_element_key"))


def _gene_metadata(mdata) -> pd.DataFrame:
    mapping = _build_gene_name_map(mdata)
    if mapping.empty:
        return pd.DataFrame(
            columns=["gene_id", "gene_name", "_gene_id_unversioned"]
        )
    metadata = pd.DataFrame(
        {
            "gene_id": mapping.index.astype(str),
            "gene_name": mapping.astype("string").fillna("").astype(str).to_numpy(),
        }
    ).drop_duplicates("gene_id", keep="first")
    metadata["_gene_id_unversioned"] = metadata["gene_id"].str.split(".").str[0]
    return metadata


def _rename_columns(schema_names: set[str]) -> dict[str, str]:
    aliases = {
        "log2_fc": "perturbo_log2_fc",
        "p_value": "perturbo_p_value",
        "log2_fc_std": "perturbo_fc_se",
        "q_value": "perturbo_q_value",
    }
    return {
        source: destination
        for source, destination in aliases.items()
        if source in schema_names and destination not in schema_names
    }


def _add_gene_names(frame, genes, *, fill_missing: bool):
    import polars as pl

    direct = genes.select(
        "gene_id", pl.col("gene_name").alias("_gene_name_direct")
    )
    unversioned = (
        genes.select(
            "_gene_id_unversioned",
            pl.col("gene_name").alias("_gene_name_unversioned"),
        )
        .unique("_gene_id_unversioned", keep="first", maintain_order=True)
    )
    gene_name = pl.coalesce("_gene_name_direct", "_gene_name_unversioned")
    if fill_missing:
        gene_name = gene_name.fill_null("")
    return (
        frame.with_columns(
            pl.col("gene_id")
            .cast(pl.String)
            .str.split_exact(".", 1)
            .struct.field("field_0")
            .alias("_gene_id_unversioned")
        )
        .join(direct.lazy(), on="gene_id", how="left", maintain_order="left")
        .join(
            unversioned.lazy(),
            on="_gene_id_unversioned",
            how="left",
            maintain_order="left",
        )
        .with_columns(gene_name.alias("gene_name"))
        .drop("_gene_name_direct", "_gene_name_unversioned", "_gene_id_unversioned")
    )


def _prepare_scan(input_path: str | Path):
    import polars as pl

    scan = pl.scan_parquet(input_path)
    schema = lazy_schema(scan)
    schema_names = set(schema.names())
    scan = scan.rename(_rename_columns(schema_names))
    renamed_schema = lazy_schema(scan)
    required = {"gene_id", "perturbo_log2_fc", "perturbo_p_value"}
    missing = required.difference(renamed_schema.names())
    if missing:
        raise ValueError(
            f"Fast Parquet enrichment is missing required columns: {sorted(missing)}"
        )
    if "perturbo_q_value" not in renamed_schema.names():
        raise ValueError(
            "Fast Parquet enrichment requires precomputed perturbo_q_value. "
            "Use the pandas fallback when q-values are absent."
        )
    if "perturbo_fdr_log10_p_value" in renamed_schema.names():
        scan = scan.drop("perturbo_fdr_log10_p_value")
        renamed_schema = lazy_schema(scan)
    if "perturbo_fc_se" not in renamed_schema.names():
        scan = scan.with_columns(
            pl.lit(None, dtype=pl.Float64).alias("perturbo_fc_se")
        )

    # Current PerTurbo Parquet outputs contain complete q-values. Checking the
    # single q-value column is cheap compared with constructing 211M enriched
    # rows and protects the existing fill-missing-q semantics.
    has_missing_q = (
        scan.select(pl.col("perturbo_q_value").is_null().any())
        .collect(engine="streaming")
        .item()
    )
    if has_missing_q:
        raise ValueError(
            "Fast Parquet enrichment does not fill missing q-values; use the "
            "pandas fallback for this input."
        )

    p_dtype = lazy_schema(scan)["perturbo_p_value"]
    scan = scan.with_columns(
        (
            pl.col("perturbo_p_value")
            .cast(pl.Float64)
            .clip(lower_bound=P_VALUE_FLOOR)
            .log10()
            .neg()
            .cast(p_dtype)
            .alias("perturbo_negLog10p")
        )
    )
    return scan


def enrich_global_parquet(
    input_path: str | Path,
    output_path: str | Path,
    mdata,
    table_kind: Literal["guide", "element"],
) -> Path:
    """Enrich a global PerTurbo table with bounded-memory Polars execution."""

    import polars as pl

    input_path = Path(input_path)
    output_path = Path(output_path)
    scan = _prepare_scan(input_path).with_columns(pl.col("gene_id").cast(pl.String))
    genes = pl.from_pandas(_gene_metadata(mdata))

    if table_kind == "guide":
        scan = scan.with_columns(pl.col("guide_id").cast(pl.String))
        metadata = pl.from_pandas(_guide_metadata(mdata))
        scan = scan.join(
            metadata.lazy(), on="guide_id", how="left", maintain_order="left"
        )
        scan = _add_gene_names(scan, genes, fill_missing=True)
        output_columns = GUIDE_OUTPUT_COLUMNS
    elif table_kind == "element":
        scan = scan.with_columns(
            pl.col("intended_target_name").cast(pl.String).str.strip_chars(),
            pl.col("intended_target_chr").cast(pl.String).str.strip_chars(),
            pl.col("intended_target_start").cast(pl.Int64, strict=False),
            pl.col("intended_target_end").cast(pl.Int64, strict=False),
        )
        metadata = pl.from_pandas(_element_metadata(mdata))
        scan = scan.join(
            metadata.lazy(),
            on=ELEMENT_COLUMNS,
            how="left",
            maintain_order="left",
        )
        scan = scan.with_columns(
            pl.col("intended_target_name").alias("element_id"),
            pl.col("intended_target_name").alias("element_name"),
            pl.col("intended_target_chr").alias("element_chr"),
            pl.col("intended_target_start").alias("element_start"),
            pl.col("intended_target_end").alias("element_end"),
        )
        scan = _add_gene_names(scan, genes, fill_missing=False)
        output_columns = ELEMENT_OUTPUT_COLUMNS
    else:
        raise ValueError(f"Unknown table kind: {table_kind}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"Streaming {table_kind} enrichment with Polars "
        f"({pl.thread_pool_size()} threads): {input_path} -> {output_path}",
        flush=True,
    )
    scan.select(output_columns).sink_parquet(
        output_path,
        compression="zstd",
        maintain_order=True,
        engine="streaming",
    )
    print(f"Completed streaming {table_kind} enrichment: {output_path}", flush=True)
    return output_path


def supports_fast_parquet_path(*paths: str | Path) -> bool:
    if not all(str(path).lower().endswith((".parquet", ".pq")) for path in paths):
        return False
    try:
        import polars  # noqa: F401
        import pyarrow  # noqa: F401
    except ImportError:
        return False
    return True
