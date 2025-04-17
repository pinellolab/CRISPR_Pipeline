#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import muon as mu
import numpy as np 

def plot_umi_threshold(mudata, save_dir):
    plt.figure(figsize=(10, 6))
    plt.rcParams['font.size'] = 12

    all_cutoffs = [200, 500, 1000, 2000, 5000]
    cell_numbers = [(mudata.mod['gene'].X.sum(1) > i).sum() for i in all_cutoffs]

    df_threshold = pd.DataFrame({
        'cell_number': cell_numbers,
        'limit': all_cutoffs
    })

    ax = sns.barplot(x='limit', y='cell_number', data=df_threshold, palette="deep")

    for i, v in enumerate(df_threshold['cell_number']):
        ax.text(i, v, str(v), ha='center', va='bottom')

    plt.title('Number of scRNA barcodes using different\nTotal UMI thresholds', fontsize=14, fontweight='bold')
    plt.xlabel('Total UMI threshold', fontsize=12)
    plt.ylabel('Number of cells', fontsize=12)

    plt.xticks(rotation=45, ha='right')
    ax.set_xticklabels(all_cutoffs)
    plt.ylim(0, max(cell_numbers) * 1.1)
    sns.despine()

    plt.tight_layout()
    plot_path = os.path.join(save_dir, 'scRNA_barcodes_UMI_thresholds.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def plot_detected_genes_threshold(mudata, save_dir):
    plt.figure(figsize=(10, 6))
    plt.rcParams['font.size'] = 12

    all_cutoffs = [200, 500, 1000, 1500, 2000]
    cell_numbers = [((mudata.mod['gene'].X > 0).sum(1) > i).sum() for i in all_cutoffs]

    df_threshold = pd.DataFrame({
        'cell_number': cell_numbers,
        'limit': all_cutoffs
    })

    ax = sns.barplot(x='limit', y='cell_number', data=df_threshold, palette="deep")

    for i, v in enumerate(df_threshold['cell_number']):
        ax.text(i, v, str(v), ha='center', va='bottom')

    plt.title('Number of scRNA barcodes with the number of detected genes > threshold', fontsize=14, fontweight='bold')
    plt.xlabel('Total detected genes threshold', fontsize=12)
    plt.ylabel('Number of cells', fontsize=12)

    plt.xticks(rotation=45, ha='right')
    ax.set_xticklabels(all_cutoffs)
    plt.ylim(0, max(cell_numbers) * 1.1)
    sns.despine()

    plt.tight_layout()
    plot_path = os.path.join(save_dir, 'scRNA_barcodes_detected_genes_thresholds.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def plot_guides_umi_threshold(mudata, save_dir):
    plt.figure(figsize=(10, 6))
    plt.rcParams['font.size'] = 12

    all_cutoffs = [0, 5, 10, 50, 100, 1000]
    cell_numbers = [(mudata.mod['guide'].X.sum(1) > i).sum() for i in all_cutoffs]

    df_threshold = pd.DataFrame({
        'cell_number': cell_numbers,
        'limit': all_cutoffs
    })

    ax = sns.barplot(x='limit', y='cell_number', data=df_threshold, palette="deep")

    for i, v in enumerate(df_threshold['cell_number']):
        ax.text(i, v, str(v), ha='center', va='bottom')

    plt.title('Number of Cells Using Different Total Guides UMI Thresholds', fontsize=14, fontweight='bold')
    plt.xlabel('Guides Total UMI Threshold', fontsize=12)
    plt.ylabel('Number of Cells', fontsize=12)

    plt.xticks(rotation=45, ha='right')
    ax.set_xticklabels(all_cutoffs)
    plt.ylim(0, max(cell_numbers) * 1.1)
    sns.despine()

    plt.tight_layout()
    plot_path = os.path.join(save_dir, 'guides_UMI_thresholds.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def plot_sgRNA_frequencies(mudata, save_dir):
    cell_ids = mudata.mod['guide'].obs.index
    guide_ids = mudata.mod['guide'].var.index
    guide_assignment_matrix = mudata.mod['guide'].layers['guide_assignment'].toarray()
    df_guide_assignment = pd.DataFrame(guide_assignment_matrix, index=cell_ids, columns=guide_ids)

    plt.figure(figsize=(10, 6))
    plt.rcParams['font.size'] = 12

    sgRNA_frequencies = (df_guide_assignment > 0).sum(axis=0)
    df_sgRNA_frequencies = sgRNA_frequencies.reset_index()
    df_sgRNA_frequencies.columns = ['sgRNA', 'Frequency']
    df_sgRNA_frequencies = df_sgRNA_frequencies.sort_values(by='Frequency', ascending=True)

    colors = sns.color_palette("Spectral", len(df_sgRNA_frequencies))
    sns.barplot(x='sgRNA', y='Frequency', data=df_sgRNA_frequencies, palette=colors)

    plt.title('Histogram of the number of sgRNA represented at the guide barcodes', fontsize=14, fontweight='bold')
    plt.xlabel('sgRNAs Represented on the Guide Barcode \n (unfiltered)', fontsize=12)
    plt.ylabel('Number of Guide Barcodes with a Given sgRNA', fontsize=12)

    plt.xticks(rotation=90, fontsize=4)
    sns.despine()

    plt.tight_layout()
    plot_path = os.path.join(save_dir, 'guides_hist_num_sgRNA.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def plot_guides_per_cell(mudata, save_dir):
    guides_per_cell = np.sum(mudata.mod['guide'].X, axis=1)
    plt.figure(figsize=(10, 6))
    plt.hist(guides_per_cell, bins=50, alpha=0.6, color='skyblue', ec="steelblue")
    plt.xlabel('Number of Guides per Cell')
    plt.ylabel('Density')
    plt.title('Histogram of Guides per Cell')
    plot_path = os.path.join(save_dir, 'guides_per_cell_histogram.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def plot_cells_per_guide(mudata, save_dir):
    cells_per_guide = np.sum(mudata.mod['guide'].X, axis=0)
    sorted_cells_per_guide = np.sort(np.array(cells_per_guide).flatten())
    plt.figure(figsize=(10, 6))
    plt.hist(sorted_cells_per_guide, bins=30, alpha=0.6, color='skyblue', ec='steelblue')
    plt.xlabel('Number of Cells per Guide')
    plt.ylabel('Density')
    plt.title('Histogram of Cells per Guide')
    plot_path = os.path.join(save_dir, 'cells_per_guide_histogram.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Generate various plots from MuData")
    parser.add_argument('--mudata', required=True, help='Path to the mudata object')
    parser.add_argument('--output_dir', required=True, help='Directory where plots will be saved')

    args = parser.parse_args()

    # Load mudata
    mudata = mu.read_h5mu(args.mudata)
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate all plots
    plot_umi_threshold(mudata, args.output_dir)
    plot_detected_genes_threshold(mudata, args.output_dir)
    plot_guides_umi_threshold(mudata, args.output_dir)
    plot_sgRNA_frequencies(mudata, args.output_dir)
    plot_guides_per_cell(mudata, args.output_dir)
    plot_cells_per_guide(mudata, args.output_dir)

if __name__ == "__main__":
    main()
