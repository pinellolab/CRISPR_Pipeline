#!/usr/bin/env python

import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def main(adata_rna, gname_rna, min_genes, min_cells, pct_mito, reference):

    gname_rna_path = os.path.join(gname_rna, 'counts_unfiltered/cells_x_genes.genes.names.txt')
    gene_df = pd.read_csv(gname_rna_path, header=None)
    gene_names = gene_df[0].tolist()

    adata_rna = sc.read(adata_rna)

    if len(gene_names) == adata_rna.shape[1]:
        # keep ensembl id as var_names
        # matched gene names as column: gene_id
        # adata_rna.var_names = gene_names
        adata_rna.var["symbol"] = gene_names
        # modify the ensembl id in adata_rna.var_names
        adata_rna.var_names = adata_rna.var_names.str.split('.').str[0]
    else:
        raise ValueError("The number of gene names does not match the number of variables in adata_rna")

    # Save knee plot
    if not os.path.exists('figures'):
        os.makedirs('figures')
        print(f"Directory '{'figures'}' created.")
    else:
        print(f"Directory already exists.")

    # knee plots
    knee_df = pd.DataFrame({
        'sum': np.array(adata_rna.X.sum(1)).flatten(),
        'barcodes': adata_rna.obs_names.values})
    knee_df = knee_df.sort_values('sum', ascending=False).reset_index(drop=True)
    knee_df['sum_log'] = np.log1p(knee_df['sum'])

    plt.figure(figsize=(8, 5))
    plt.plot(knee_df['sum_log'], knee_df.index, marker='o', linestyle='-', markersize=3)
    plt.ylabel('Barcode number')
    plt.xlabel('Log of UMI Counts')
    plt.title('Knee Plot')
    plt.savefig('figures/knee_plot_scRNA.png')

    # Add batch number
    adata_rna.obs['batch_number'] = adata_rna.obs['batch'].factorize()[0] + 1

    mt_prefix = "MT-" if reference == "human" else "Mt-"
    adata_rna.var["mt"] = adata_rna.var['symbol'].str.startswith('MT-')
    adata_rna.var["ribo"] = adata_rna.var['symbol'].str.startswith(("RPS", "RPL"))

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata_rna, qc_vars=["mt", "ribo"], inplace=True, log1p=True)

    # Plot violin
    sc.pl.violin(
        adata_rna,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        save='plot_scrna.png'
    )

    # Plot scatter
    sc.pl.scatter(
        adata_rna,
        x="total_counts",
        y="n_genes_by_counts",
        color="pct_counts_mt",
        save='plot_scrna.png'
    )

    # filter for mic_cells and min_genes
    sc.pp.filter_cells(adata_rna, min_genes=min_genes)
    sc.pp.filter_genes(adata_rna, min_cells=min_cells)

    # filter for percent mito
    pct_mito=pct_mito 
    adata_rna = adata_rna[adata_rna.obs['pct_counts_mt'] < pct_mito, :]

    # Save
    adata_rna.write('filtered_anndata.h5ad')

    print("Quality control completed and saved to 'filtered_anndata.h5ad'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC on AnnData.')
    parser.add_argument('adata_rna', type=str, help='Path to the AnnData file.')
    parser.add_argument('gname_rna', type=str, help='Path to the cells x genes txt file.')
    parser.add_argument('--min_genes', type=int, default=100, help='Minimum number of genes per cell.')
    parser.add_argument('--min_cells', type=int, default=3, help='Minimum number of cells per gene.')
    parser.add_argument('--pct_mito', type=float, default=0.2, help='Minimum percent of proportion of mitochondrial reads in cells.')
    parser.add_argument('--reference', type=str, required=True, help='Reference species')

    args = parser.parse_args()
    main(args.adata_rna, args.gname_rna, args.min_genes, args.min_cells, args.pct_mito, args.reference)
