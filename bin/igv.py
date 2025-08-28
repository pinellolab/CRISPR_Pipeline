#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from collections import defaultdict
import muon as mu
import os
from gtfparse import read_gtf
from typing import Literal, Tuple, Dict, Optional

def process_coordinates(mdata) -> Dict[str, list]:
    """Extract coordinate information from MuData object"""
    coordinate_dict = {}

    # Process gene coordinates
    for index, row in mdata.mod["gene"].var.iterrows():
        if np.isnan(row["gene_start"]) or np.isnan(row["gene_end"]):
            continue
        coordinate_dict[index] = [row["gene_chr"], row["gene_start"], row["gene_end"]]

    # Process guide coordinates
    for index, row in mdata.mod["guide"].var.iterrows():
        if row["intended_target_name"] in coordinate_dict or row["intended_target_name"] == "non-targeting":
            continue
        coordinate_dict[row["intended_target_name"]] = [
            row["intended_target_chr"],
            row["intended_target_start"],
            row["intended_target_end"]
        ]

    return coordinate_dict

def igv(mdata, gtf: str, method: Optional[Literal['sceptre', 'perturbo']] = None,
        results_key: str = 'test_results') -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Generate bedpe and bedgraph data for generic or method-specific data"""
    # Set column names based on whether method is specified
    if not method:  # Generic columns
        log2_fc_col = "log2_fc"
        p_value_col = "p_value"
    else:  # Method-specific columns
        log2_fc_col = f"{method}_log2_fc"
        p_value_col = f"{method}_p_value"

    # Process coordinates
    coordinate_dict = process_coordinates(mdata)

    # Process GTF file
    df_gtf = read_gtf(gtf).to_pandas()
    gencode_df = df_gtf[['gene_id', 'gene_name']]
    gencode_df['gene_id2'] = gencode_df['gene_id'].str.split('.').str[0]
    gencode_df = gencode_df.drop_duplicates()

    # Initialize data structures
    bedpe = defaultdict(list)
    bedgraph = defaultdict(list)

    # Process test results
    test_results = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})
    merged_df = test_results.merge(
        gencode_df[['gene_id2', 'gene_name']],
        left_on='gene_id',
        right_on='gene_id2',
        how='left'
    )

    # Filter out rows where required columns are missing
    merged_df = merged_df.dropna(subset=[log2_fc_col, p_value_col])

    for index, row in merged_df.iterrows():
        if row["intended_target_name"] == row["gene_name"]:
            # PROMOTER interactions
            if row["intended_target_name"] in coordinate_dict:
                coords = coordinate_dict[row["intended_target_name"]]
                bedgraph["chr"].append(coords[0])
                bedgraph["start"].append(coords[1])
                bedgraph["end"].append(coords[2])
                bedgraph["p_value"].append(row[p_value_col])
                bedgraph["log2_fc"].append(row[log2_fc_col])
        else:
            # ENHANCER-GENE interactions
            if (row["intended_target_name"] in coordinate_dict and
                row["gene_id"] in coordinate_dict):
                source_coords = coordinate_dict[row["intended_target_name"]]
                target_coords = coordinate_dict[row["gene_id"]]

                bedpe["chr1"].append(source_coords[0])
                bedpe["start1"].append(source_coords[1])
                bedpe["end1"].append(source_coords[2])
                bedpe["chr2"].append(target_coords[0])
                bedpe["start2"].append(target_coords[1])
                bedpe["end2"].append(target_coords[2])
                bedpe["p_value"].append(row[p_value_col])
                bedpe["log2_fc"].append(row[log2_fc_col])

    bedpe_df = pd.DataFrame(bedpe)
    bedgraph_df = pd.DataFrame(bedgraph)

    method_name = method if method else "Analysis"
    if bedpe_df.empty:
        print(f"Warning: {method_name} bedpe_df is empty.")
    if bedgraph_df.empty:
        print(f"Warning: {method_name} bedgraph_df is empty.")

    print(f"\n{method_name.capitalize()} statistics ({results_key}):")
    print(f"Number of enhancer-gene interactions: {len(bedpe_df)}")
    print(f"Number of promoter interactions: {len(bedgraph_df)}")

    return bedpe_df, bedgraph_df

def process_results_config(mdata, gtf: str, results_key: str, analysis_type: Optional[str] = None):
    """Process a single results configuration (either test_results or cis/trans_test_results)"""

    # Check if the key exists in mdata.uns
    if results_key not in mdata.uns:
        print(f"Warning: {results_key} not found in mdata.uns, skipping...")
        return

    print(f"\nProcessing {results_key}...")

    # Check available methods/columns
    results_df = pd.DataFrame(mdata.uns[results_key])
    cols = results_df.columns
    print(f"Available columns for {results_key}:", cols)

    # Create output directory
    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)

    # Check for generic columns first
    if 'log2_fc' in cols and 'p_value' in cols:
        print(f"Using generic log2_fc and p_value columns for {results_key}")
        bedpe_df, bedgraph_df = igv(mdata, gtf, None, results_key)

        # Save generic files with analysis type prefix
        filename_prefix = f"{analysis_type}_" if analysis_type else ""
        bedpe_path = os.path.join(output_dir, f"{filename_prefix}evaluation.bedpe")
        bedgraph_path = os.path.join(output_dir, f"{filename_prefix}evaluation.bedgraph")

        bedpe_df.to_csv(bedpe_path, sep="\t", index=False, header=False)
        bedgraph_df.to_csv(bedgraph_path, sep="\t", index=False, header=False)

        print(f"\nFiles saved for {results_key}:")
        print(f"bedpe file: {bedpe_path}")
        print(f"bedgraph file: {bedgraph_path}")
    else:
        # Check for method-specific columns
        available_methods = []
        if 'sceptre_log2_fc' in cols and not results_df['sceptre_log2_fc'].isna().all():
            available_methods.append('sceptre')
        if 'perturbo_log2_fc' in cols and not results_df['perturbo_log2_fc'].isna().all():
            available_methods.append('perturbo')

        print(f"Available methods for {results_key}: {available_methods}")

        if not available_methods:
            print(f"No methods with valid data found for {results_key}")
            return

        # Process data for available methods
        for method in available_methods:
            bedpe_df, bedgraph_df = igv(mdata, gtf, method, results_key)

            # Generate outputs with analysis type prefix
            filename_prefix = f"{analysis_type}_" if analysis_type else ""
            bedpe_path = os.path.join(output_dir, f"{filename_prefix}{method}.bedpe")
            bedgraph_path = os.path.join(output_dir, f"{filename_prefix}{method}.bedgraph")

            # Save files
            bedpe_df.to_csv(bedpe_path, sep="\t", index=False, header=False)
            bedgraph_df.to_csv(bedgraph_path, sep="\t", index=False, header=False)

            print(f"\n{method.capitalize()} files saved for {results_key}:")
            print(f"bedpe file: {bedpe_path}")
            print(f"bedgraph file: {bedgraph_path}")

if __name__ == "__main__":
    print("Starting program...")
    parser = argparse.ArgumentParser(description="Process MuData and generate bedpe and bedgraph files")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--gtf", type=str, required=True, help="Path to the GTF file")
    parser.add_argument("--results_key", type=str, default="test_results",
                      help="Key for test results in mdata.uns")
    parser.add_argument("--default", action="store_true",
                      help="Process mudata with cis_per_element_results and trans_per_element_results instead of single test_results")

    args = parser.parse_args()

    print("Loading MuData file...")
    mdata = mu.read(args.mdata_path)

    # Determine which results to process based on --default flag
    if args.default:
        # Process both cis and trans results
        results_configs = [
            {"key": "cis_per_element_results", "type": "cis"},
            {"key": "trans_per_element_results", "type": "trans"}
        ]
    else:
        # Process single test_results
        results_configs = [
            {"key": args.results_key, "type": None}
        ]

    # Process each results configuration
    for config in results_configs:
        process_results_config(mdata, args.gtf, config["key"], config["type"])

    print("\nAll processing completed successfully")
