#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np
from scipy.stats import false_discovery_control
from analysis_output_formatting import (
    format_element_output,
    format_guide_output,
    make_h5mu_safe_dataframe,
)


ELEMENT_BASE_KEYS = ["gene_id", "intended_target_name"]
ELEMENT_LOCATION_KEYS = [
    "intended_target_chr",
    "intended_target_start",
    "intended_target_end",
]


def _build_merge_keys(sceptre_df: pd.DataFrame, perturbo_df: pd.DataFrame):
    full_keys = ELEMENT_BASE_KEYS + ELEMENT_LOCATION_KEYS
    if all(col in sceptre_df.columns for col in full_keys) and all(
        col in perturbo_df.columns for col in full_keys
    ):
        return full_keys, True
    return ELEMENT_BASE_KEYS, False


def _with_merge_key_columns(df: pd.DataFrame, key_cols):
    out = df.copy()
    merge_cols = []
    for col in key_cols:
        merge_col = f"__merge_{col}"
        out[merge_col] = out[col].astype("string").fillna("__NA__")
        merge_cols.append(merge_col)
    return out, merge_cols


def _bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if not valid.any():
        return out

    out.loc[p.loc[valid].index] = false_discovery_control(
        p.loc[valid].to_numpy(dtype=float), method="bh"
    )
    return out


def _add_perturbo_q_value(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "perturbo_p_value" not in out.columns:
        return out

    out["perturbo_q_value"] = _bh_adjust(out["perturbo_p_value"])
    return out


def _add_perturbo_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(
        columns={
            "log2_fc": "perturbo_log2_fc",
            "p_value": "perturbo_p_value",
            "log2_fc_std": "perturbo_fc_se",
            "q_value": "perturbo_q_value",
        }
    )
    out = out.drop(columns=["perturbo_fdr_log10_p_value"], errors="ignore")

    if "perturbo_p_value" in out.columns:
        out["perturbo_q_value"] = _bh_adjust(out["perturbo_p_value"])
    if "perturbo_fc_se" not in out.columns:
        out["perturbo_fc_se"] = np.nan

    return out


def _add_sceptre_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.rename(
        columns={
            "log2_fc": "sceptre_log2_fc",
            "p_value": "sceptre_p_value",
            "q_value": "sceptre_q_value",
            "se_fold_change": "sceptre_fc_se",
        }
    )

    if "sceptre_p_value" in out.columns:
        computed_sceptre_q = _bh_adjust(out["sceptre_p_value"])
        if "sceptre_q_value" not in out.columns:
            out["sceptre_q_value"] = computed_sceptre_q
        else:
            out["sceptre_q_value"] = out["sceptre_q_value"].fillna(
                computed_sceptre_q
            )

    if "sceptre_fc_se" not in out.columns:
        out["sceptre_fc_se"] = np.nan

    return out


def _existing_columns(df: pd.DataFrame, columns):
    return [col for col in columns if col in df.columns]


def merge_method_results(sceptre_per_guide, sceptre_per_element, perturbo_per_guide, perturbo_per_element, base_mudata_path):
    """
    Merge SCEPTRE and PerTurbo results into a single MuData object.
    
    Args:
        sceptre_per_guide: Path to SCEPTRE per_guide_output.tsv
        sceptre_per_element: Path to SCEPTRE per_element_output.tsv  
        perturbo_per_guide: Path to PerTurbo per_guide_output.tsv
        perturbo_per_element: Path to PerTurbo per_element_output.tsv
        base_mudata_path: Path to base mudata file for structure
    """
    print("Loading input files...")
    
    # Load SCEPTRE results
    sceptre_guide_df = pd.read_csv(sceptre_per_guide, sep='\t')
    sceptre_element_df = pd.read_csv(sceptre_per_element, sep='\t')
    
    # Rename SCEPTRE columns to indicate method
    sceptre_guide_df = _add_sceptre_columns(sceptre_guide_df)
    sceptre_element_df = _add_sceptre_columns(sceptre_element_df)
    
    # Load PerTurbo results
    perturbo_guide_df = pd.read_csv(perturbo_per_guide, sep='\t')
    perturbo_element_df = pd.read_csv(perturbo_per_element, sep='\t')
    
    # Rename PerTurbo columns to indicate method
    perturbo_guide_df = _add_perturbo_columns(perturbo_guide_df)
    perturbo_element_df = _add_perturbo_columns(perturbo_element_df)
    
    print("Merging per-guide results...")
    # Merge per-guide results
    sceptre_guide_cols = _existing_columns(
        sceptre_guide_df,
        [
            "gene_id",
            "guide_id",
            "sceptre_log2_fc",
            "sceptre_p_value",
            "sceptre_q_value",
            "sceptre_fc_se",
        ],
    )
    perturbo_guide_cols = _existing_columns(
        perturbo_guide_df,
        [
            "gene_id",
            "guide_id",
            "perturbo_log2_fc",
            "perturbo_p_value",
            "perturbo_q_value",
            "perturbo_fc_se",
        ],
    )
    merged_guide_df = pd.merge(
        sceptre_guide_df[sceptre_guide_cols],
        perturbo_guide_df[perturbo_guide_cols],
        on=['gene_id', 'guide_id'],
        how='outer'
    )
    merged_guide_df = _add_perturbo_q_value(merged_guide_df)
    
    print("Merging per-element results...")
    merge_keys, using_full_coordinate_keys = _build_merge_keys(
        sceptre_element_df,
        perturbo_element_df,
    )
    if using_full_coordinate_keys:
        print("Using coordinate-aware element merge keys:", ", ".join(merge_keys))
    else:
        print(
            "Falling back to legacy per-element merge keys: gene_id, intended_target_name"
        )

    sceptre_element_cols = _existing_columns(
        sceptre_element_df,
        merge_keys
        + [
            "sceptre_log2_fc",
            "sceptre_p_value",
            "sceptre_q_value",
            "sceptre_fc_se",
        ],
    )
    perturbo_element_cols = merge_keys + [
        "perturbo_log2_fc",
        "perturbo_p_value",
        "perturbo_q_value",
        "perturbo_fc_se",
    ]
    sceptre_element_merge_df = sceptre_element_df[sceptre_element_cols].copy()
    perturbo_element_merge_df = perturbo_element_df[perturbo_element_cols].copy()

    sceptre_element_merge_df, left_merge_cols = _with_merge_key_columns(
        sceptre_element_merge_df,
        merge_keys,
    )
    perturbo_element_merge_df, right_merge_cols = _with_merge_key_columns(
        perturbo_element_merge_df,
        merge_keys,
    )

    merged_element_df = pd.merge(
        sceptre_element_merge_df,
        perturbo_element_merge_df,
        left_on=left_merge_cols,
        right_on=right_merge_cols,
        how='outer',
        suffixes=('', '_perturbo'),
    )

    for key_col in merge_keys:
        right_key_col = f"{key_col}_perturbo"
        if right_key_col in merged_element_df.columns:
            merged_element_df[key_col] = merged_element_df[key_col].combine_first(
                merged_element_df[right_key_col]
            )
            merged_element_df.drop(columns=[right_key_col], inplace=True)

    merged_element_df.drop(columns=left_merge_cols + right_merge_cols, inplace=True, errors='ignore')

    merged_element_df = _add_perturbo_q_value(merged_element_df)

    preferred_order = _existing_columns(
        merged_element_df,
        merge_keys
        + [
            "sceptre_log2_fc",
            "sceptre_p_value",
            "sceptre_q_value",
            "sceptre_fc_se",
            "perturbo_log2_fc",
            "perturbo_p_value",
            "perturbo_q_value",
            "perturbo_fc_se",
        ],
    )
    merged_element_df = merged_element_df[preferred_order]
    
    # Load base mudata for structure
    base_mdata = mu.read_h5mu(base_mudata_path)

    merged_guide_df = format_guide_output(merged_guide_df)
    merged_element_df = format_element_output(merged_element_df, base_mdata)
    
    # Store merged results in mudata
    base_mdata.uns['per_guide_results'] = make_h5mu_safe_dataframe(merged_guide_df)
    base_mdata.uns['per_element_results'] = make_h5mu_safe_dataframe(merged_element_df)
    
    # Write outputs
    print("Writing merged results...")
    base_mdata.write("inference_mudata.h5mu", compression="gzip")
    merged_guide_df.to_csv("per_guide_output.tsv.gz", sep='\t', index=False, compression='gzip')
    merged_element_df.to_csv("per_element_output.tsv.gz", sep='\t', index=False, compression='gzip')
    
    print("Successfully merged results from both methods!")
    return base_mdata

def main():
    parser = argparse.ArgumentParser(description='Merge SCEPTRE and PerTurbo results')
    parser.add_argument('--sceptre_per_guide', required=True, help='Path to SCEPTRE per_guide_output.tsv')
    parser.add_argument('--sceptre_per_element', required=True, help='Path to SCEPTRE per_element_output.tsv')
    parser.add_argument('--perturbo_per_guide', required=True, help='Path to PerTurbo per_guide_output.tsv')
    parser.add_argument('--perturbo_per_element', required=True, help='Path to PerTurbo per_element_output.tsv')
    parser.add_argument('--base_mudata', required=True, help='Path to base mudata file for structure')
    
    args = parser.parse_args()
    
    merge_method_results(
        args.sceptre_per_guide,
        args.sceptre_per_element, 
        args.perturbo_per_guide,
        args.perturbo_per_element,
        args.base_mudata
    )

if __name__ == "__main__":
    main()
