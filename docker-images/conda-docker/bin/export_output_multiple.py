#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np
from scipy import sparse
import time

def merge_results(sceptre_result, perturbo_mudata):
    """
    Merge SCEPTRE and Perturbo results and save to MuData.
    """
    start_time = time.time()
    print("Loading and merging results...")
    
    # Load and rename SCEPTRE results
    test_result = pd.read_csv(sceptre_result)
    test_result.rename(columns={'log2_fc': 'sceptre_log2_fc', 'p_value': 'sceptre_p_value'}, inplace=True)

    # Load and rename Perturbo results
    mudata = mu.read_h5mu(perturbo_mudata)
    perturbo_result = pd.DataFrame(mudata.uns['test_results'])
    
    # Merge results efficiently (only include necessary columns from perturbo_result)
    merged_result = pd.merge(
        test_result,
        perturbo_result[['guide_id', 'gene_id', 'perturbo_log2_fc', 'perturbo_p_value']],
        on=['guide_id', 'gene_id'],
        how='left'  # Use left join to keep all SCEPTRE results
    )

    # Update MuData with merged results
    mudata.uns["test_results"] = merged_result.to_dict(orient='list')
    mudata.write("inference_mudata.h5mu")
    print(f"Updated MuData saved to 'inference_mudata.h5mu' in {time.time() - start_time:.2f} seconds.")

    return merged_result, mudata

def export_output(merged_result, mudata):
    """
    Generate and save per-guide and per-element outputs.
    Optimized for performance.
    """
    start_time = time.time()
    print("Generating outputs...")
    
    # Generate per-guide output more efficiently
    guide_var = mudata.mod['guide'].var.reset_index()
    per_guide_output = pd.merge(
        merged_result,
        guide_var,
        how='left', 
        on=['intended_target_name','guide_id']
    ).rename(columns={'guide_id': 'guide_id(s)'})
    
    # OPTIMIZATION: Vectorized gene expression calculation
    print("Calculating average gene expressions...")
    
    # Create mapping from gene symbol to matrix column indices
    gene_X = mudata.mod['gene'].X
    gene_symbols = mudata.mod['gene'].var['symbol'].values
    
    # Create a dictionary to map gene symbols to their column indices
    symbol_to_indices = {}
    for idx, symbol in enumerate(gene_symbols):
        if symbol in symbol_to_indices:
            symbol_to_indices[symbol].append(idx)
        else:
            symbol_to_indices[symbol] = [idx]
    
    # Calculate average expression for each unique target gene
    unique_targets = per_guide_output['intended_target_name'].unique()
    target_to_avg_exp = {}
    
    # Check if we have a sparse matrix
    is_sparse = sparse.issparse(gene_X)
    
    for target in unique_targets:
        if target in symbol_to_indices and symbol_to_indices[target]:
            indices = symbol_to_indices[target]
            if is_sparse:
                # For sparse matrices, extract only the required columns
                subset = gene_X[:, indices]
                avg_exp = np.mean(subset.toarray())
            else:
                # For dense matrices, use direct indexing
                avg_exp = np.mean(gene_X[:, indices])
        else:
            avg_exp = np.nan
        target_to_avg_exp[target] = avg_exp
    
    # Apply the precomputed averages
    per_guide_output['avg_gene_expression'] = per_guide_output['intended_target_name'].map(target_to_avg_exp)
    
    # Add cell number column
    per_guide_output['cell_number'] = mudata.shape[0]

    # Reorder columns
    per_guide_output = per_guide_output.reindex(columns=[
        'intended_target_name', 'guide_id(s)', 'intended_target_chr',
        'intended_target_start', 'intended_target_end', 'gene_id',
        'sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc',
        'perturbo_p_value', 'cell_number', 'avg_gene_expression'
    ])
    
    # OPTIMIZATION: More efficient per-element output generation
    print("Generating per-element output...")
    
    # First, create a mapping of gene_id to concatenated guide_ids
    guide_id_map = per_guide_output.groupby('gene_id')['guide_id(s)'].apply(
        lambda x: ','.join(x.dropna().astype(str))
    ).to_dict()
    
    # Then, get one row per gene_id as a template
    per_element_template = per_guide_output.drop_duplicates('gene_id').drop(columns=['guide_id(s)'])
    
    # Create per-element dataframe efficiently
    per_element_output = pd.DataFrame({'gene_id': list(guide_id_map.keys())})
    per_element_output['guide_id(s)'] = per_element_output['gene_id'].map(guide_id_map)
    
    # Merge with template
    per_element_output = pd.merge(
        per_element_output,
        per_element_template,
        on='gene_id',
        how='left'
    )

    # Reorder columns
    column_order = [
        'intended_target_name',
        'guide_id(s)',
        'intended_target_chr',
        'intended_target_start',
        'intended_target_end',
        'gene_id',
        'sceptre_log2_fc',
        'sceptre_p_value',
        'perturbo_log2_fc',
        'perturbo_p_value',
        'cell_number',
        'avg_gene_expression'
    ]
    per_element_output = per_element_output[column_order]

    # Save outputs
    print("Saving output files...")
    per_guide_output.to_csv('per_guide_output.tsv', sep='\t', index=False)
    per_element_output.to_csv('per_element_output.tsv', sep='\t', index=False)
    
    print(f"Exported output files successfully in {time.time() - start_time:.2f} seconds.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process and merge SCEPTRE and Perturbo results.')
    parser.add_argument('--sceptre_result', required=True, help='Path to the SCEPTRE results CSV file.')
    parser.add_argument('--perturbo_mudata', required=True, help='Path to the Perturbo h5mu file.')
    args = parser.parse_args()

    total_start_time = time.time()
    # Process results and export outputs
    merged_result, mudata = merge_results(args.sceptre_result, args.perturbo_mudata)
    export_output(merged_result, mudata)
    print(f"Total processing time: {time.time() - total_start_time:.2f} seconds.")