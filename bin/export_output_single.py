#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
import numpy as np
from scipy import sparse
import time

def export_output(mudata, inference_method):
    """
    Generate per-guide and per-element outputs from MuData and save as TSV files.
    Optimized for performance.
    """
    start_time = time.time()
    
    # Check if test_results exists in uns
    if 'test_results' not in mudata.uns:
        print("Warning: No test results found in mudata.uns")
        return

    # Get test results 
    test_results = pd.DataFrame(mudata.uns['test_results'])
    if test_results.empty:
        print("Warning: Test results dataframe is empty")
        return

    # Define column mappings
    col_map = {
        'sceptre': {'log2_fc': 'sceptre_log2_fc', 'p_value': 'sceptre_p_value'},
        'perturbo': {'log2_fc': 'perturbo_log2_fc', 'p_value': 'perturbo_p_value'}
    }
    
    if inference_method not in col_map:
        raise ValueError("Invalid inference_method. Must be 'sceptre' or 'perturbo'.")

    print("Processing data...")
    try:
        # Merge with guide data 
        guide_var_df = mudata.mod['guide'].var.reset_index()
        
        per_guide_output = pd.merge(
            test_results,
            guide_var_df,
            how='left',
            on=['intended_target_name', 'guide_id']
        )

        # Pre-extract gene expression matrix and symbol array for faster lookup
        print("Preparing gene expression data...")
        gene_X = mudata.mod['gene'].X
        gene_symbols = mudata.mod['gene'].var['symbol'].values
        
        # Convert to array once if sparse 
        is_sparse = sparse.issparse(gene_X)
        if is_sparse:
            print("Gene expression matrix is sparse, optimizing calculations...")
        
        # Create a mapping from gene symbol to column index for fast lookup
        symbol_to_idx = {}
        for idx, symbol in enumerate(gene_symbols):
            if symbol in symbol_to_idx:
                symbol_to_idx[symbol].append(idx)
            else:
                symbol_to_idx[symbol] = [idx]
        
        # OPTIMIZATION: Calculate average gene expression vectorized where possible
        print("Calculating average gene expression...")
        unique_targets = per_guide_output['intended_target_name'].unique()
        target_to_avg_exp = {}
        
        for target in unique_targets:
            if target in symbol_to_idx:
                col_indices = symbol_to_idx[target]
                if len(col_indices) > 0:
                    if is_sparse:
                        # Extract only the columns we need
                        subset = gene_X[:, col_indices]
                        avg_exp = subset.mean()
                    else:
                        avg_exp = np.mean(gene_X[:, col_indices])
                else:
                    avg_exp = np.nan
            else:
                avg_exp = np.nan
            target_to_avg_exp[target] = avg_exp
        
        # Map the averages back to the dataframe
        per_guide_output['avg_gene_expression'] = per_guide_output['intended_target_name'].map(target_to_avg_exp)

        # Add cell number column
        per_guide_output['cell_number'] = mudata.shape[0]

        # Rename columns based on inference method
        per_guide_output = per_guide_output.rename(columns=col_map[inference_method])

        # OPTIMIZATION: Add missing inference method columns all at once
        missing_cols = {
            'sceptre_log2_fc': None, 
            'sceptre_p_value': None, 
            'perturbo_log2_fc': None, 
            'perturbo_p_value': None
        }
        for col, default in missing_cols.items():
            if col not in per_guide_output.columns:
                per_guide_output[col] = default

        # Reorder and rename columns
        final_cols = [
            'intended_target_name', 'guide_id', 'intended_target_chr',
            'intended_target_start', 'intended_target_end', 'gene_id',
            'sceptre_log2_fc', 'sceptre_p_value', 'perturbo_log2_fc',
            'perturbo_p_value', 'cell_number', 'avg_gene_expression'
        ]

        per_guide_output = per_guide_output[final_cols].rename(columns={'guide_id': 'guide_id(s)'})

        # Generate per-element output
        print("Generating per-element output...")
        
        # OPTIMIZATION: Faster groupby approach
        # First, create a mapping of gene_id to concatenated guide_ids
        guide_id_map = per_guide_output.groupby('gene_id')['guide_id(s)'].apply(
            lambda x: ','.join(x.dropna().astype(str))
        ).to_dict()
        
        # Then, get one row per gene_id as a template
        per_element_template = per_guide_output.drop_duplicates('gene_id').drop(columns=['guide_id(s)'])
        
        # Create per-element dataframe
        per_element_output = pd.DataFrame({'gene_id': list(guide_id_map.keys())})
        per_element_output['guide_id(s)'] = per_element_output['gene_id'].map(guide_id_map)
        
        # Merge with template
        per_element_output = pd.merge(
            per_element_output,
            per_element_template,
            on='gene_id',
            how='left'
        )

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
        
        end_time = time.time()
        print(f"Exported output files successfully in {end_time - start_time:.2f} seconds.")
        
    except Exception as e:
        print(f"Error during export: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Export outputs from MuData.")
    parser.add_argument("--inference_method", type=str, required=True, help="Inference method: 'sceptre' or 'perturbo'.")
    parser.add_argument("--mudata", type=str, required=True, help="Path to input MuData file.")
    args = parser.parse_args()

    try:
        print(f"Loading MuData file from {args.mudata}...")
        start_time = time.time()
        mudata = mu.read_h5mu(args.mudata)
        load_time = time.time() - start_time
        print(f"MuData loaded in {load_time:.2f} seconds.")
        
        export_output(mudata, args.inference_method)
    except Exception as e:
        print(f"Error loading or processing MuData file: {str(e)}")
        raise