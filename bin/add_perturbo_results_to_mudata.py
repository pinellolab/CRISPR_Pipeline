#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np
from analysis_output_formatting import make_h5mu_safe_dataframe
from result_table_io import read_result_table

def add_perturbo_results_to_mudata(per_guide_results, per_element_results, base_mudata_path, output_path):
    """
    Add PerTurbo per_guide and per_element results to a base MuData file.
    
    Args:
        per_guide_results: Path to per-guide TSV or Parquet results
        per_element_results: Path to per-element TSV or Parquet results
        base_mudata_path: Path to base mudata file
        output_path: Output path for updated MuData file
    """
    print("Loading input files...")
    
    per_guide_df = read_result_table(per_guide_results)
    per_element_df = read_result_table(per_element_results)
    
    # Load base mudata
    mdata = mu.read_h5mu(base_mudata_path)
    
    print("Adding results to MuData...")
    
    # Store results in mudata .uns field
    mdata.uns['per_guide_results'] = make_h5mu_safe_dataframe(per_guide_df)
    mdata.uns['per_element_results'] = make_h5mu_safe_dataframe(per_element_df)
    
    # Write the updated mudata
    print(f"Writing updated MuData to {output_path}...")
    mdata.write(output_path, compression="gzip")
    
    print("Successfully added PerTurbo results to MuData!")
    print(f"Output contains:")
    print(f"  - per_guide_results: {len(per_guide_df)} entries")
    print(f"  - per_element_results: {len(per_element_df)} entries")
    
    return mdata

def main():
    parser = argparse.ArgumentParser(description='Add PerTurbo results to MuData file')
    parser.add_argument(
        '--per_guide_results',
        '--per_guide_tsv',
        dest='per_guide_results',
        required=True,
        help='Path to per-guide .tsv[.gz] or .parquet results',
    )
    parser.add_argument(
        '--per_element_results',
        '--per_element_tsv',
        dest='per_element_results',
        required=True,
        help='Path to per-element .tsv[.gz] or .parquet results',
    )
    parser.add_argument('--base_mudata', required=True, help='Path to base mudata file')
    parser.add_argument('--output', required=True, help='Output path for updated MuData file')
    
    args = parser.parse_args()
    
    add_perturbo_results_to_mudata(
        args.per_guide_results,
        args.per_element_results,
        args.base_mudata,
        args.output
    )

if __name__ == "__main__":
    main()
