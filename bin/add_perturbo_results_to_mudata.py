#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np

def add_perturbo_results_to_mudata(per_guide_tsv, per_element_tsv, base_mudata_path, output_path):
    """
    Add PerTurbo per_guide and per_element results to a base MuData file.
    
    Args:
        per_guide_tsv: Path to per_guide_output.tsv
        per_element_tsv: Path to per_element_output.tsv
        base_mudata_path: Path to base mudata file
        output_path: Output path for updated MuData file
    """
    print("Loading input files...")
    
    # Load TSV files
    per_guide_df = pd.read_csv(per_guide_tsv, sep='\t')
    per_element_df = pd.read_csv(per_element_tsv, sep='\t')
    
    # Load base mudata
    mdata = mu.read_h5mu(base_mudata_path)
    
    print("Adding results to MuData...")
    
    # Store results in mudata .uns field
    mdata.uns['per_guide_results'] = per_guide_df
    mdata.uns['per_element_results'] = per_element_df
    
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
    parser.add_argument('--per_guide_tsv', required=True, help='Path to per_guide_output.tsv')
    parser.add_argument('--per_element_tsv', required=True, help='Path to per_element_output.tsv')
    parser.add_argument('--base_mudata', required=True, help='Path to base mudata file')
    parser.add_argument('--output', required=True, help='Output path for updated MuData file')
    
    args = parser.parse_args()
    
    add_perturbo_results_to_mudata(
        args.per_guide_tsv,
        args.per_element_tsv,
        args.base_mudata,
        args.output
    )

if __name__ == "__main__":
    main()
