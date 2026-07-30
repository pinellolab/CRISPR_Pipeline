#!/usr/bin/env python3

import argparse
from typing import Iterable, List, Optional

import mudata as mu
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control


JOIN_COLUMNS = ["gene_id", "guide_id"]
GUIDE_METADATA_COLUMNS = [
    "guide_sequence",
    "guide_type",
    "targeting",
    "guide_chr",
    "guide_start",
    "guide_end",
    "guide_strand",
    "pam",
    "intended_target_name",
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
]
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
    "guide_id",
] + GUIDE_METADATA_COLUMNS + ["gene_name", "gene_id", "nPerturbedCells"]
P_VALUE_FLOOR = 1e-300


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
    return -np.log10(values.where(values.notna()).clip(lower=pvalue_floor))


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if valid.any():
        out.loc[valid] = false_discovery_control(
            p.loc[valid].to_numpy(dtype=float), method="bh"
        )
    return out


def _assert_unique_keys(df: pd.DataFrame, keys: List[str], label: str) -> None:
    duplicates = df.duplicated(keys, keep=False)
    if duplicates.any():
        examples = df.loc[duplicates, keys].head(10).to_dict(orient="records")
        raise ValueError(
            f"{label} contains duplicate rows for key columns {keys}. "
            f"Sample duplicates: {examples}"
        )


def _build_guide_metadata(guide) -> pd.DataFrame:
    guide_var = guide.var.copy()
    if "guide_id" not in guide_var.columns:
        guide_var["guide_id"] = guide_var.index.astype(str)
    guide_var["guide_id"] = guide_var["guide_id"].astype(str)
    _assert_unique_keys(guide_var, ["guide_id"], "guide metadata")

    alias_map = {
        "guide_sequence": ["guide_sequence", "spacer", "sequence"],
        "guide_type": ["guide_type", "type"],
        "targeting": ["targeting"],
        "guide_chr": ["guide_chr"],
        "guide_start": ["guide_start"],
        "guide_end": ["guide_end"],
        "guide_strand": ["guide_strand", "strand"],
        "pam": ["pam"],
        "intended_target_name": ["intended_target_name"],
        "intended_target_chr": ["intended_target_chr"],
        "intended_target_start": ["intended_target_start"],
        "intended_target_end": ["intended_target_end"],
    }
    metadata = pd.DataFrame({"guide_id": guide_var["guide_id"]})
    for output_col, candidates in alias_map.items():
        source_col = _first_existing_column_optional(guide_var, candidates)
        metadata[output_col] = (
            guide_var[source_col].to_numpy() if source_col is not None else pd.NA
        )

    assignment = (
        guide.layers["guide_assignment"]
        if "guide_assignment" in guide.layers
        else guide.X
    )
    if sparse.issparse(assignment):
        counts = np.asarray((assignment > 0).sum(axis=0)).ravel()
    else:
        counts = np.asarray(assignment > 0).sum(axis=0).ravel()
    if len(counts) != len(metadata):
        raise ValueError(
            "Guide assignment columns do not match the number of guide metadata rows."
        )
    metadata["nPerturbedCells"] = pd.array(counts, dtype="Int64")
    return metadata


def _build_gene_name_map(mdata) -> pd.Series:
    for modality_name in ("gene", "rna"):
        if not hasattr(mdata, "mod") or modality_name not in mdata.mod:
            continue
        gene_var = mdata[modality_name].var.copy()
        if gene_var.empty:
            continue
        name_col = (
            "symbol"
            if "symbol" in gene_var.columns
            else "gene_name"
            if "gene_name" in gene_var.columns
            else None
        )
        if name_col is None:
            continue
        mapping = pd.Series(
            gene_var[name_col].astype("string").values,
            index=gene_var.index.astype(str),
        )
        return mapping[~mapping.index.duplicated(keep="first")]
    return pd.Series(dtype="string")


def _fill_gene_names(results: pd.DataFrame, mdata) -> pd.Series:
    mapping = _build_gene_name_map(mdata)
    if mapping.empty:
        return pd.Series(pd.NA, index=results.index, dtype="string")
    out = results["gene_id"].astype(str).map(mapping)
    stripped_map = pd.Series(
        mapping.values, index=mapping.index.to_series().str.split(".").str[0]
    )
    stripped_map = stripped_map[~stripped_map.index.duplicated(keep="first")]
    missing = out.isna()
    out.loc[missing] = (
        results.loc[missing, "gene_id"].astype(str).str.split(".").str[0].map(stripped_map)
    )
    return out


def _prepare_metric_subset(
    results: pd.DataFrame,
    method: str,
    fc_candidates: List[str],
    p_candidates: List[str],
    q_candidates: List[str],
    se_candidates: List[str],
    label: str,
) -> pd.DataFrame:
    normalized = results.copy()
    normalized["gene_id"] = normalized["gene_id"].astype(str)
    normalized["guide_id"] = normalized["guide_id"].astype(str)
    fc_col = _first_existing_column(normalized, fc_candidates, f"{label} log2_fc")
    p_col = _first_existing_column(normalized, p_candidates, f"{label} p_value")
    q_col = _first_existing_column_optional(normalized, q_candidates)
    se_col = _first_existing_column_optional(normalized, se_candidates)

    metric_cols = [fc_col, p_col] + [c for c in (q_col, se_col) if c is not None]
    subset = normalized[JOIN_COLUMNS + metric_cols].copy()
    _assert_unique_keys(subset, JOIN_COLUMNS, label)
    subset = subset.rename(
        columns={
            fc_col: f"{method}_log2_fc",
            p_col: f"_{method}_p_value",
            **({q_col: f"{method}_q_value"} if q_col else {}),
            **({se_col: f"{method}_fc_se"} if se_col else {}),
        }
    )
    computed_q = _bh_adjust(subset[f"_{method}_p_value"])
    q_output = f"{method}_q_value"
    subset[q_output] = (
        subset[q_output].fillna(computed_q)
        if q_output in subset.columns
        else computed_q
    )
    se_output = f"{method}_fc_se"
    if se_output not in subset.columns:
        subset[se_output] = np.nan
    return subset


def create_catalog_per_guide(
    local_analysis_per_guide: pd.DataFrame,
    global_analysis_per_guide: pd.DataFrame,
    mdata,
    pvalue_floor: float = P_VALUE_FLOOR,
) -> pd.DataFrame:
    local = _prepare_metric_subset(
        local_analysis_per_guide,
        "sceptre",
        ["sceptre_log2_fc", "log2_fc"],
        ["sceptre_p_value", "p_value"],
        ["sceptre_q_value", "q_value"],
        ["sceptre_fc_se", "se_fold_change"],
        "local_analysis_per_guide",
    )
    global_results = _prepare_metric_subset(
        global_analysis_per_guide,
        "perturbo",
        ["perturbo_log2_fc", "log2_fc"],
        ["perturbo_p_value", "p_value"],
        ["perturbo_q_value", "q_value"],
        ["perturbo_fc_se", "log2_fc_std"],
        "global_analysis_per_guide",
    )
    merged = local.merge(global_results, on=JOIN_COLUMNS, how="outer")
    merged = merged.merge(_build_guide_metadata(mdata["guide"]), on="guide_id", how="left")

    for method in ("sceptre", "perturbo"):
        raw = merged[f"_{method}_p_value"]
        merged[f"{method}_p_value"] = raw
        merged[f"{method}_negLog10p"] = _neg_log10(raw, pvalue_floor)
    merged["gene_name"] = _fill_gene_names(merged, mdata)

    catalog = merged[OUTPUT_COLUMNS].copy()
    return catalog.sort_values(
        ["guide_chr", "guide_start", "guide_end", "guide_id", "gene_id"],
        kind="stable",
        na_position="last",
    ).reset_index(drop=True)


def build_catalog_per_guide_output(
    local_analysis_per_guide_path: str,
    global_analysis_per_guide_path: str,
    mudata_path: str,
    output_path: str,
    pvalue_floor: float = P_VALUE_FLOOR,
) -> pd.DataFrame:
    local_results = pd.read_csv(local_analysis_per_guide_path, sep="\t")
    global_results = pd.read_csv(global_analysis_per_guide_path, sep="\t")
    mdata = mu.read_h5mu(mudata_path)
    catalog = create_catalog_per_guide(
        local_results, global_results, mdata, pvalue_floor=pvalue_floor
    )
    catalog.to_csv(output_path, sep="\t", index=False, compression="gzip")
    return catalog


def main():
    parser = argparse.ArgumentParser(
        description="Build the catalog per-guide table from local and global analysis results."
    )
    parser.add_argument("--local_analysis_per_guide", required=True)
    parser.add_argument("--global_analysis_per_guide", required=True)
    parser.add_argument("--mudata", required=True)
    parser.add_argument("--output", default="catalog_per_guide_output.tsv.gz")
    parser.add_argument("--pvalue_floor", type=float, default=P_VALUE_FLOOR)
    args = parser.parse_args()
    build_catalog_per_guide_output(
        local_analysis_per_guide_path=args.local_analysis_per_guide,
        global_analysis_per_guide_path=args.global_analysis_per_guide,
        mudata_path=args.mudata,
        output_path=args.output,
        pvalue_floor=args.pvalue_floor,
    )


if __name__ == "__main__":
    main()
