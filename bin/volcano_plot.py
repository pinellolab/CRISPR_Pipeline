#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
import muon as mu
from mudata import MuData
import matplotlib.pyplot as plt
import os
from typing import Literal, Optional

def volcano(mdata: MuData, method: Optional[Literal['sceptre', 'perturbo']],
           log2_fc_thresh: float, p_value_thresh: float, ax: plt.Axes = None,
           results_key: str = 'test_results', analysis_type: Optional[str] = None):
    if ax is None:
        ax = plt.gca()

    # Get test results from MuData
    test_results = pd.DataFrame({k: v for k, v in mdata.uns[results_key].items()})

    # Determine which columns to use
    if not method:  # If method is None, try to use generic columns
        if 'log2_fc' in test_results.columns and 'p_value' in test_results.columns:
            log2_fc_col = 'log2_fc'
            p_value_col = 'p_value'
        else:
            raise ValueError("No generic log2_fc and p_value columns found")
    else:  # Use method-specific columns
        log2_fc_col = f"{method}_log2_fc"
        p_value_col = f"{method}_p_value"

    # Plot all data points (non-significant genes)
    ax.scatter(x=test_results[log2_fc_col],
              y=test_results[p_value_col].apply(lambda x: -np.log10(x)),
              s=5, color="green", label="Not significant")

    # Highlight down or up-regulated genes based on the thresholds
    down = test_results[(test_results[log2_fc_col] <= -log2_fc_thresh) &
                       (test_results[p_value_col] <= p_value_thresh)]
    up = test_results[(test_results[log2_fc_col] >= log2_fc_thresh) &
                      (test_results[p_value_col] <= p_value_thresh)]

    filtered_down = down[np.isfinite(down[log2_fc_col]) &
                        np.isfinite(down[p_value_col])]

    method_name = method if method else "Analysis"

    down_copy = down.copy()
    down_copy[log2_fc_col] = down_copy[log2_fc_col].replace([np.inf, -np.inf], [1e10, -1e10])
    down_copy[p_value_col] = down_copy[p_value_col].replace([np.inf, -np.inf], [1e10, -1e10])

    top_5_down_log2fc = down_copy.sort_values(by=log2_fc_col, ascending=True).head(10)
    top_5_down_pval = down_copy.sort_values(by=p_value_col, ascending=True).head(10)

    annotated = pd.concat([top_5_down_log2fc, top_5_down_pval])

    # Plot down-regulated genes
    ax.scatter(x=down[log2_fc_col],
              y=down[p_value_col].apply(lambda x: -np.log10(x)),
              s=10, label="Down-regulated", color="blue")

    # Plot up-regulated genes
    ax.scatter(x=up[log2_fc_col],
              y=up[p_value_col].apply(lambda x: -np.log10(x)),
              s=10, label="Up-regulated", color="red")

    # Annotate genes with offset labels
    for index, row in annotated.iterrows():
        x = row[log2_fc_col]
        y = -np.log10(row[p_value_col])
        label = f"{row['intended_target_name']} ({row['guide_id']})"

        # Add small offset to avoid overlapping with points
        offset_x = 0.1 if x >= 0 else -0.1
        offset_y = 0.1

        # Add annotation with a simple arrow
        ax.annotate(label,
                   xy=(x, y),
                   xytext=(x + offset_x, y + offset_y),
                   fontsize=8,
                   arrowprops=dict(arrowstyle='->', color='gray', lw=0.5))

    # Set plot boundaries and labels
    low, high = ax.get_xlim()
    bound = max(abs(low), abs(high)) + 0.5

    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("p-value (-log10)")

    # Draw threshold lines
    ax.axvline(-log2_fc_thresh, color="grey", linestyle="--")
    ax.axvline(log2_fc_thresh, color="grey", linestyle="--")
    ax.axhline(-np.log10(p_value_thresh), color="grey", linestyle="--")

    ax.set_xlim(-bound, bound)
    ax.legend(loc="upper right")

    method_title = method.capitalize() if method else "Differential Expression"
    title_prefix = f"{analysis_type}_" if analysis_type else ""
    ax.set_title(f"{title_prefix}{method_title} Volcano Plot")

def process_results_config(mdata: MuData, args, results_key: str, analysis_type: Optional[str] = None):
    """Process a single results configuration (either test_results or cis/trans_test_results)"""

    # Check if the key exists in mdata.uns
    if results_key not in mdata.uns:
        print(f"Warning: {results_key} not found in mdata.uns, skipping...")
        return

    print(f"\nProcessing {results_key}...")

    # Check available methods
    results_df = pd.DataFrame(mdata.uns[results_key])
    cols = results_df.columns

    output_dir = "evaluation_output"
    os.makedirs(output_dir, exist_ok=True)

    # Check for generic columns first
    if 'log2_fc' in cols and 'p_value' in cols:
        print(f"Using generic log2_fc and p_value columns for {results_key}")
        fig, ax = plt.subplots(1, 1, figsize=(10, 10))
        volcano(mdata, None, args.log2_fc, args.p_value, ax, results_key, analysis_type)

        # Add analysis type to figure title
        fig_title = f"{analysis_type.capitalize()} Differential Expression Results" if analysis_type else "Differential Expression Results"
        plt.suptitle(fig_title, y=1.02)

        # Save plot with analysis type prefix
        filename_prefix = f"{analysis_type}_" if analysis_type else ""
        output_file = os.path.join(output_dir, f"{filename_prefix}volcano_plot.png")
        print(f"Saving plot to {output_file}")
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

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

        if len(available_methods) == 1:
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            volcano(mdata, available_methods[0], args.log2_fc, args.p_value, ax, results_key, analysis_type)

            # Add analysis type to figure title
            method_title = available_methods[0].title()
            fig_title = f"{analysis_type.capitalize()} {method_title} Results" if analysis_type else f"{method_title} Results"
            plt.suptitle(fig_title, y=1.02)

            # Save single plot with analysis type prefix
            filename_prefix = f"{analysis_type}_" if analysis_type else ""
            output_file = os.path.join(output_dir, f"{filename_prefix}{available_methods[0].lower()}_volcano_plot.png")
            print(f"Saving plot to {output_file}")
            plt.tight_layout()
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()

        elif len(available_methods) == 2:
            # Create and save Sceptre plot
            fig1, ax1 = plt.subplots(1, 1, figsize=(10, 10))
            volcano(mdata, 'sceptre', args.log2_fc, args.p_value, ax1, results_key, analysis_type)

            # Add analysis type to figure title
            fig_title = f"{analysis_type.capitalize()} Sceptre Results" if analysis_type else "Sceptre Results"
            plt.suptitle(fig_title, y=1.02)

            filename_prefix = f"{analysis_type}_" if analysis_type else ""
            sceptre_output = os.path.join(output_dir, f"{filename_prefix}sceptre_volcano_plot.png")
            print(f"Saving Sceptre plot to {sceptre_output}")
            plt.tight_layout()
            plt.savefig(sceptre_output, dpi=300, bbox_inches='tight')
            plt.close()

            # Create and save Perturbo plot
            fig2, ax2 = plt.subplots(1, 1, figsize=(10, 10))
            volcano(mdata, 'perturbo', args.log2_fc, args.p_value, ax2, results_key, analysis_type)

            # Add analysis type to figure title
            fig_title = f"{analysis_type.capitalize()} Perturbo Results" if analysis_type else "Perturbo Results"
            plt.suptitle(fig_title, y=1.02)

            perturbo_output = os.path.join(output_dir, f"{filename_prefix}perturbo_volcano_plot.png")
            print(f"Saving Perturbo plot to {perturbo_output}")
            plt.tight_layout()
            plt.savefig(perturbo_output, dpi=300, bbox_inches='tight')
            plt.close()

        print(f"Plot(s) for {results_key} saved successfully")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate comparison volcano plots from MuData")
    parser.add_argument("mdata_path", type=str, help="Path to the MuData file")
    parser.add_argument("--log2_fc", type=float, default=1,
                      help="log2 fold change threshold")
    parser.add_argument("--p_value", type=float, default=0.01,
                      help="p-value threshold")
    parser.add_argument("--results_key", type=str, default="test_results",
                      help="Key for test results in mdata.uns")
    parser.add_argument("--default", action="store_true",
                      help="Process mudata with cis_test_results and trans_test_results instead of single test_results")

    args = parser.parse_args()

    mdata = mu.read(args.mdata_path)

    # Determine which results to process based on --default flag
    if args.default:
        # Process both cis and trans results
        results_configs = [
            {"key": "cis_test_results", "type": "cis"},
            {"key": "trans_test_results", "type": "trans"}
        ]
    else:
        # Process single test_results
        results_configs = [
            {"key": args.results_key, "type": None}
        ]

    # Process each results configuration
    for config in results_configs:
        process_results_config(mdata, args, config["key"], config["type"])

    print("All plots generated successfully")
