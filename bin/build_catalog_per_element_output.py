#!/usr/bin/env python3

import argparse
from typing import Dict, Iterable, List, Optional

import mudata as mu
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from scipy import sparse
from result_table_io import categoricalize_text_columns, read_result_table, write_result_table

ELEMENT_COLUMNS = [
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
]
JOIN_COLUMNS = ["gene_id"] + ELEMENT_COLUMNS
OUTPUT_COLUMNS = [
    "sceptre_log2_fc",
    "sceptre_p_value",
    "sceptre_q_value",
    "sceptre_fc_se",
    "sceptre_negLog10p",
    "perturbo_log2_fc",
    "perturbo_p_value",
    "perturbo_q_value",
    "perturbo_fc_se",
    "perturbo_negLog10p",
    "element_id",
    "element_type",
    "element_chr",
    "element_start",
    "element_end",
    "element_name",
    "guide_ids",
    "num_guides",
    "gene_name",
    "gene_id",
    "nPerturbedCells",
]
P_VALUE_FLOOR = 1e-300
MISSING_TOKEN = ""


def _normalize_text(series: pd.Series) -> pd.Series:
    return series.astype("string").str.strip()


def _normalize_int(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def _normalize_element_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in ELEMENT_COLUMNS:
        if col not in out.columns:
            out[col] = pd.NA

    out["intended_target_name"] = _normalize_text(out["intended_target_name"])
    out["intended_target_chr"] = _normalize_text(out["intended_target_chr"])
    out["intended_target_start"] = _normalize_int(out["intended_target_start"])
    out["intended_target_end"] = _normalize_int(out["intended_target_end"])
    return out


def _element_key(df: pd.DataFrame) -> pd.Series:
    key_parts = [
        _normalize_text(df["intended_target_name"]).fillna(MISSING_TOKEN),
        _normalize_text(df["intended_target_chr"]).fillna(MISSING_TOKEN),
        df["intended_target_start"].astype("string").fillna(MISSING_TOKEN),
        df["intended_target_end"].astype("string").fillna(MISSING_TOKEN),
    ]
    return key_parts[0] + "|" + key_parts[1] + "|" + key_parts[2] + "|" + key_parts[3]


def _first_existing_column(
    df: pd.DataFrame, candidates: Iterable[str], label: str
) -> str:
    for col in candidates:
        if col in df.columns:
            return col
    raise KeyError(f"Missing {label} column. Tried: {list(candidates)}")


def _first_existing_column_optional(
    df: pd.DataFrame, candidates: Iterable[str]
) -> Optional[str]:
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _neg_log10(series: pd.Series, pvalue_floor: float) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    finite = values.where(values.notna())
    clamped = finite.clip(lower=pvalue_floor)
    return -np.log10(clamped)


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    """
    Benjamini-Hochberg FDR correction.
    NaN inputs remain NaN.
    """
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if not valid.any():
        return out

    q = false_discovery_control(p.loc[valid].to_numpy(dtype=float), method="bh")
    out.loc[p.loc[valid].index] = q
    return out


def _assert_unique_keys(df: pd.DataFrame, keys: List[str], label: str) -> None:
    dup_mask = df.duplicated(keys, keep=False)
    if dup_mask.any():
        examples = df.loc[dup_mask, keys].head(10).to_dict(orient="records")
        raise ValueError(
            f"{label} contains duplicate rows for key columns {keys}. "
            f"Sample duplicates: {examples}"
        )


def _collapse_unique_strings(values: pd.Series) -> str:
    cleaned = sorted(
        {str(v).strip() for v in values if pd.notna(v) and str(v).strip() != ""}
    )
    if not cleaned:
        return pd.NA
    return ";".join(cleaned)


def _build_element_metadata(guide) -> pd.DataFrame:
    guide_var = guide.var.copy()
    if "guide_id" not in guide_var.columns:
        guide_var["guide_id"] = guide_var.index.astype(str)

    guide_var = _normalize_element_columns(guide_var)
    if "type" not in guide_var.columns:
        guide_var["type"] = pd.NA

    guide_var["guide_id"] = guide_var["guide_id"].astype("string")
    guide_var["_element_key"] = _element_key(guide_var)

    grouped = guide_var.groupby(
        ELEMENT_COLUMNS + ["_element_key"], dropna=False, as_index=False
    ).agg(
        guide_ids=("guide_id", _collapse_unique_strings),
        num_guides=("guide_id", "nunique"),
        element_type=("type", _collapse_unique_strings),
    )
    grouped["num_guides"] = grouped["num_guides"].astype("Int64")

    assignment = (
        guide.layers["guide_assignment"]
        if "guide_assignment" in guide.layers
        else guide.X
    )
    if not sparse.issparse(assignment):
        assignment = sparse.csr_matrix(assignment)
    else:
        assignment = assignment.tocsr()

    element_codes, unique_keys = pd.factorize(guide_var["_element_key"], sort=False)
    n_guides = len(element_codes)
    indicator = sparse.csr_matrix(
        (
            np.ones(n_guides, dtype=np.int8),
            (np.arange(n_guides), element_codes),
        ),
        shape=(n_guides, len(unique_keys)),
    )

    cell_by_element = assignment @ indicator
    n_perturbed = np.asarray((cell_by_element > 0).sum(axis=0)).ravel().astype(int)
    n_perturbed_map: Dict[str, int] = {
        key: int(count)
        for key, count in zip(unique_keys.tolist(), n_perturbed.tolist())
    }

    grouped["nPerturbedCells"] = (
        grouped["_element_key"].map(n_perturbed_map).astype("Int64")
    )
    return grouped


def _build_gene_name_map(mdata) -> pd.Series:
    gene_var = mdata["gene"].var.copy()
    if gene_var is None or gene_var.empty:
        return pd.Series(dtype="string")

    if "symbol" in gene_var.columns:
        name_col = "symbol"
    elif "gene_name" in gene_var.columns:
        name_col = "gene_name"
    else:
        return pd.Series(dtype="string")

    mapping = pd.Series(
        gene_var[name_col].astype("string").values, index=gene_var.index.astype(str)
    )
    mapping = mapping[~mapping.index.duplicated(keep="first")]
    return mapping


def _fill_gene_names(results: pd.DataFrame, mdata) -> pd.Series:
    gene_name_map = _build_gene_name_map(mdata)
    if gene_name_map.empty:
        return pd.Series(pd.NA, index=results.index, dtype="string")

    out = results["gene_id"].astype(str).map(gene_name_map)
    if out.notna().all():
        return out

    stripped_index = gene_name_map.index.to_series().str.split(".").str[0]
    stripped_map = (
        pd.DataFrame(
            {"gene_id": stripped_index.values, "gene_name": gene_name_map.values}
        )
        .dropna(subset=["gene_id"])
        .drop_duplicates(subset=["gene_id"], keep="first")
        .set_index("gene_id")["gene_name"]
    )

    missing_mask = out.isna()
    out.loc[missing_mask] = (
        results.loc[missing_mask, "gene_id"]
        .astype(str)
        .str.split(".")
        .str[0]
        .map(stripped_map)
    )
    return out


def create_catalog_per_element(
    local_analysis_per_element: pd.DataFrame,
    global_analysis_per_element: pd.DataFrame,
    mdata,
    pvalue_floor: float = P_VALUE_FLOOR,
) -> pd.DataFrame:
    local_results = _normalize_element_columns(local_analysis_per_element)
    global_results = _normalize_element_columns(global_analysis_per_element)

    local_results["gene_id"] = local_results["gene_id"].astype(str)
    global_results["gene_id"] = global_results["gene_id"].astype(str)

    cis_fc_col = _first_existing_column(
        local_results, ["sceptre_log2_fc", "log2_fc"], "local-analysis SCEPTRE log2_fc"
    )
    cis_p_col = _first_existing_column(
        local_results, ["sceptre_p_value", "p_value"], "local-analysis SCEPTRE p_value"
    )
    cis_q_col = _first_existing_column_optional(
        local_results, ["sceptre_q_value", "q_value"]
    )
    cis_fc_se_col = _first_existing_column_optional(
        local_results, ["sceptre_fc_se", "se_fold_change"]
    )
    trans_fc_col = _first_existing_column(
        global_results, ["perturbo_log2_fc", "log2_fc"], "global-analysis PerTurbo log2_fc"
    )
    trans_p_col = _first_existing_column(
        global_results, ["perturbo_p_value", "p_value"], "global-analysis PerTurbo p_value"
    )
    trans_q_col = _first_existing_column_optional(
        global_results, ["perturbo_q_value", "q_value"]
    )
    trans_fc_se_col = _first_existing_column_optional(
        global_results, ["perturbo_fc_se", "log2_fc_std"]
    )

    cis_metric_cols = [cis_fc_col, cis_p_col]
    if cis_q_col is not None:
        cis_metric_cols.append(cis_q_col)
    if cis_fc_se_col is not None:
        cis_metric_cols.append(cis_fc_se_col)

    cis_subset = local_results[JOIN_COLUMNS + cis_metric_cols].copy()
    trans_metric_cols = [trans_fc_col, trans_p_col]
    if trans_q_col is not None:
        trans_metric_cols.append(trans_q_col)
    if trans_fc_se_col is not None:
        trans_metric_cols.append(trans_fc_se_col)

    trans_subset = global_results[JOIN_COLUMNS + trans_metric_cols].copy()
    _assert_unique_keys(cis_subset, JOIN_COLUMNS, "local_analysis_per_element")
    _assert_unique_keys(trans_subset, JOIN_COLUMNS, "global_analysis_per_element")

    cis_rename_cols = {
        cis_fc_col: "sceptre_log2_fc",
        cis_p_col: "_sceptre_p_value",
    }
    if cis_q_col is not None:
        cis_rename_cols[cis_q_col] = "sceptre_q_value"
    if cis_fc_se_col is not None:
        cis_rename_cols[cis_fc_se_col] = "sceptre_fc_se"
    cis_subset = cis_subset.rename(columns=cis_rename_cols)
    computed_sceptre_q = _bh_adjust(cis_subset["_sceptre_p_value"])
    if "sceptre_q_value" not in cis_subset.columns:
        cis_subset["sceptre_q_value"] = computed_sceptre_q
    else:
        cis_subset["sceptre_q_value"] = cis_subset["sceptre_q_value"].fillna(
            computed_sceptre_q
        )
    if "sceptre_fc_se" not in cis_subset.columns:
        cis_subset["sceptre_fc_se"] = np.nan

    trans_rename_cols = {
        trans_fc_col: "perturbo_log2_fc",
        trans_p_col: "_perturbo_p_value",
    }
    if trans_q_col is not None:
        trans_rename_cols[trans_q_col] = "perturbo_q_value"
    if trans_fc_se_col is not None:
        trans_rename_cols[trans_fc_se_col] = "perturbo_fc_se"
    trans_subset = trans_subset.rename(columns=trans_rename_cols)
    computed_perturbo_q = _bh_adjust(trans_subset["_perturbo_p_value"])
    if "perturbo_q_value" not in trans_subset.columns:
        trans_subset["perturbo_q_value"] = computed_perturbo_q
    else:
        trans_subset["perturbo_q_value"] = trans_subset["perturbo_q_value"].fillna(
            computed_perturbo_q
        )
    if "perturbo_fc_se" not in trans_subset.columns:
        trans_subset["perturbo_fc_se"] = np.nan

    merged = cis_subset.merge(trans_subset, on=JOIN_COLUMNS, how="outer")

    merged["_element_key"] = _element_key(merged)
    element_meta = _build_element_metadata(mdata["guide"])
    merged = merged.merge(
        element_meta[
            ["_element_key", "guide_ids", "num_guides", "element_type", "nPerturbedCells"]
        ],
        on="_element_key",
        how="left",
    )

    merged["sceptre_negLog10p"] = _neg_log10(merged["_sceptre_p_value"], pvalue_floor)
    merged["sceptre_p_value"] = merged["_sceptre_p_value"]
    merged["perturbo_negLog10p"] = _neg_log10(
        merged["_perturbo_p_value"], pvalue_floor
    )
    merged["perturbo_p_value"] = merged["_perturbo_p_value"]

    merged["element_name"] = merged["intended_target_name"]
    merged["element_id"] = merged["element_name"]
    merged["element_chr"] = merged["intended_target_chr"]
    merged["element_start"] = merged["intended_target_start"]
    merged["element_end"] = merged["intended_target_end"]

    merged["gene_name"] = _fill_gene_names(merged, mdata)

    catalog = categoricalize_text_columns(merged[OUTPUT_COLUMNS])
    catalog = catalog.sort_values(
        by=["element_chr", "element_start", "element_end", "element_name", "gene_id"],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)

    return catalog


def build_catalog_per_element_output(
    local_analysis_per_element_path: str,
    global_analysis_per_element_path: str,
    mudata_path: str,
    output_path: str,
    pvalue_floor: float = P_VALUE_FLOOR,
) -> pd.DataFrame:
    local_results = read_result_table(local_analysis_per_element_path)
    global_results = read_result_table(global_analysis_per_element_path)
    # backed="r" avoids loading gene/guide .X into memory; this script only
    # reads .var and .layers["guide_assignment"], which load eagerly either way.
    mdata = mu.read_h5mu(mudata_path, backed="r")

    catalog = create_catalog_per_element(
        local_results, global_results, mdata, pvalue_floor=pvalue_floor
    )
    write_result_table(catalog, output_path)
    return catalog


def main():
    parser = argparse.ArgumentParser(
        description="Build the catalog per-element table from local and global analysis results."
    )
    parser.add_argument(
        "--local_analysis_per_element",
        required=True,
        help="Path to local-analysis per-element TSV(.gz)",
    )
    parser.add_argument(
        "--global_analysis_per_element",
        required=True,
        help="Path to global-analysis per-element TSV(.gz)",
    )
    parser.add_argument(
        "--mudata", required=True, help="Path to inference MuData (h5mu)"
    )
    parser.add_argument(
        "--output",
        default="catalog_per_element_output.tsv.gz",
        help="Output path for catalog table (.tsv[.gz] or .parquet)",
    )
    parser.add_argument(
        "--pvalue_floor",
        type=float,
        default=P_VALUE_FLOOR,
        help="Lower bound used before -log10 transform.",
    )

    args = parser.parse_args()
    build_catalog_per_element_output(
        local_analysis_per_element_path=args.local_analysis_per_element,
        global_analysis_per_element_path=args.global_analysis_per_element,
        mudata_path=args.mudata,
        output_path=args.output,
        pvalue_floor=args.pvalue_floor,
    )


if __name__ == "__main__":
    main()
