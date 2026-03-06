#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu


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
    sceptre_guide_df = sceptre_guide_df.rename(columns={
        'log2_fc': 'sceptre_log2_fc',
        'p_value': 'sceptre_p_value'
    })
    sceptre_element_df = sceptre_element_df.rename(columns={
        'log2_fc': 'sceptre_log2_fc', 
        'p_value': 'sceptre_p_value'
    })
    
    # Load PerTurbo results
    perturbo_guide_df = pd.read_csv(perturbo_per_guide, sep='\t')
    perturbo_element_df = pd.read_csv(perturbo_per_element, sep='\t')
    
    # Rename PerTurbo columns to indicate method
    perturbo_guide_df = perturbo_guide_df.rename(columns={
        'log2_fc': 'perturbo_log2_fc',
        'p_value': 'perturbo_p_value'
    })
    perturbo_element_df = perturbo_element_df.rename(columns={
        'log2_fc': 'perturbo_log2_fc',
        'p_value': 'perturbo_p_value'
    })
    
    print("Merging per-guide results...")
    # Merge per-guide results
    merged_guide_df = pd.merge(
        sceptre_guide_df[['gene_id', 'guide_id', 'sceptre_log2_fc', 'sceptre_p_value']],
        perturbo_guide_df[['gene_id', 'guide_id', 'perturbo_log2_fc', 'perturbo_p_value']],
        on=['gene_id', 'guide_id'],
        how='outer'
    )
    
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

    sceptre_element_cols = merge_keys + ['sceptre_log2_fc', 'sceptre_p_value']
    perturbo_element_cols = merge_keys + ['perturbo_log2_fc', 'perturbo_p_value']
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

    preferred_order = merge_keys + ['sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc', 'perturbo_p_value']
    merged_element_df = merged_element_df[preferred_order]
    
    # Load base mudata for structure
    base_mdata = mu.read_h5mu(base_mudata_path)
    
    # Store merged results in mudata
    base_mdata.uns['per_guide_results'] = merged_guide_df
    base_mdata.uns['per_element_results'] = merged_element_df
    
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
