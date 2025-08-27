#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np
from pathlib import Path

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
        sceptre_guide_df[['gene_id', 'guide_id', 'intended_target_name', 'sceptre_log2_fc', 'sceptre_p_value']],
        perturbo_guide_df[['gene_id', 'guide_id', 'intended_target_name', 'perturbo_log2_fc', 'perturbo_p_value']],
        on=['gene_id', 'guide_id', 'intended_target_name'],
        how='outer'
    )
    
    print("Merging per-element results...")
    # Merge per-element results  
    merged_element_df = pd.merge(
        sceptre_element_df[['gene_id', 'intended_target_name', 'sceptre_log2_fc', 'sceptre_p_value']],
        perturbo_element_df[['gene_id', 'intended_target_name', 'perturbo_log2_fc', 'perturbo_p_value']],
        on=['gene_id', 'intended_target_name'],
        how='outer'
    )
    
    # Load base mudata for structure
    base_mdata = mu.read_h5mu(base_mudata_path)
    
    # Store merged results in mudata
    base_mdata.uns['per_guide_results'] = merged_guide_df
    base_mdata.uns['per_element_results'] = merged_element_df
    
    # Write outputs
    print("Writing merged results...")
    base_mdata.write("inference_mudata.h5mu", compression="gzip")
    merged_guide_df.to_csv("per_guide_output.tsv", sep='\t', index=False)
    merged_element_df.to_csv("per_element_output.tsv", sep='\t', index=False)
    
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
