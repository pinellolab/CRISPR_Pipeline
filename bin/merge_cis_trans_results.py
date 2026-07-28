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

    if "perturbo_p_value" not in out.columns:
        return out

    out["perturbo_q_value"] = _bh_adjust(out["perturbo_p_value"])
    return out


def _finalize_perturbo_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "perturbo_fc_se" not in out.columns and "perturbo_p_value" in out.columns:
        out["perturbo_fc_se"] = np.nan
    return out


def merge_cis_trans_results(cis_per_guide_path, cis_per_element_path, trans_per_guide_path, trans_per_element_path, base_mudata_path, output_path):
    """
    Merge cis and trans results from TSV files, separating results into cis/trans specific fields.
    
    Args:
        cis_per_guide_path: Path to cis per_guide_output.tsv
        cis_per_element_path: Path to cis per_element_output.tsv
        trans_per_guide_path: Path to trans per_guide_output.tsv
        trans_per_element_path: Path to trans per_element_output.tsv
        base_mudata_path: Path to base mudata file for structure
        output_path: Output path for merged MuData
    """
    print("Loading input files...")
    
    # Load TSV files
    cis_per_guide = pd.read_csv(cis_per_guide_path, sep='\t')
    cis_per_element = pd.read_csv(cis_per_element_path, sep='\t')
    trans_per_guide = pd.read_csv(trans_per_guide_path, sep='\t')
    trans_per_element = pd.read_csv(trans_per_element_path, sep='\t')

    cis_per_guide = _finalize_perturbo_columns(_add_perturbo_columns(cis_per_guide))
    cis_per_element = _finalize_perturbo_columns(_add_perturbo_columns(cis_per_element))
    trans_per_guide = _finalize_perturbo_columns(_add_perturbo_columns(trans_per_guide))
    trans_per_element = _finalize_perturbo_columns(
        _add_perturbo_columns(trans_per_element)
    )
    
    # Load base mudata for structure
    base_mdata = mu.read_h5mu(base_mudata_path)

    cis_per_guide = format_guide_output(cis_per_guide)
    trans_per_guide = format_guide_output(trans_per_guide)
    cis_per_element = format_element_output(cis_per_element, base_mdata)
    trans_per_element = format_element_output(trans_per_element, base_mdata)
    
    # Store results in separate cis/trans fields
    print("Storing separated cis/trans results...")
    base_mdata.uns['cis_per_guide_results'] = make_h5mu_safe_dataframe(cis_per_guide)
    base_mdata.uns['cis_per_element_results'] = make_h5mu_safe_dataframe(cis_per_element)
    base_mdata.uns['trans_per_guide_results'] = make_h5mu_safe_dataframe(trans_per_guide)
    base_mdata.uns['trans_per_element_results'] = make_h5mu_safe_dataframe(trans_per_element)
    
    # Remove the original combined results to avoid confusion
    if 'per_guide_results' in base_mdata.uns:
        del base_mdata.uns['per_guide_results']
    if 'per_element_results' in base_mdata.uns:
        del base_mdata.uns['per_element_results']
    
    # Write the merged mudata
    print(f"Writing merged MuData to {output_path}...")
    base_mdata.write(output_path, compression="gzip")
    
    # Write compressed TSV files with final output names
    print("Writing compressed TSV files...")
    cis_per_guide.to_csv("cis_per_guide_output.tsv.gz", sep='\t', index=False, compression='gzip')
    cis_per_element.to_csv("cis_per_element_output.tsv.gz", sep='\t', index=False, compression='gzip')
    trans_per_guide.to_csv("trans_per_guide_output.tsv.gz", sep='\t', index=False, compression='gzip')
    trans_per_element.to_csv("trans_per_element_output.tsv.gz", sep='\t', index=False, compression='gzip')
    
    print("Successfully merged cis and trans results!")
    print(f"Output contains:")
    print(f"  - cis_per_guide_results: {len(cis_per_guide)} entries")
    print(f"  - cis_per_element_results: {len(cis_per_element)} entries")
    print(f"  - trans_per_guide_results: {len(trans_per_guide)} entries")
    print(f"  - trans_per_element_results: {len(trans_per_element)} entries")
    
    return base_mdata

def main():
    parser = argparse.ArgumentParser(description='Merge cis and trans results from TSV files')
    parser.add_argument('--cis_per_guide', required=True, help='Path to cis per_guide_output.tsv')
    parser.add_argument('--cis_per_element', required=True, help='Path to cis per_element_output.tsv')
    parser.add_argument('--trans_per_guide', required=True, help='Path to trans per_guide_output.tsv')
    parser.add_argument('--trans_per_element', required=True, help='Path to trans per_element_output.tsv')
    parser.add_argument('--base_mudata', required=True, help='Path to base mudata file for structure')
    parser.add_argument('--output', required=True, help='Output path for merged MuData file')
    
    args = parser.parse_args()
    
    merge_cis_trans_results(
        args.cis_per_guide,
        args.cis_per_element,
        args.trans_per_guide,
        args.trans_per_element,
        args.base_mudata,
        args.output
    )

if __name__ == "__main__":
    main()
