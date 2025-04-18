#!/usr/bin/env python

import argparse
import anndata as ad
import pandas as pd
import numpy as np
from muon import MuData
from gtfparse import read_gtf
import matplotlib.pyplot as plt
import os

def main(adata_rna, adata_guide, adata_hashing, guide_metadata, gtf, moi, capture_method):
    # Load the data
    guide_metadata = pd.read_csv(guide_metadata, sep='\t')
    adata_rna = ad.read_h5ad(adata_rna)
    adata_guide = ad.read_h5ad(adata_guide)
    adata_hashing = ad.read_h5ad(adata_hashing)
    df_gtf = read_gtf(gtf).to_pandas()

    ## change in adata_guide
    # adding var for guide
    adata_guide.var.reset_index(inplace=True)
    # rename gene_id in guide
    adata_guide.var.rename(columns={'gene_id': 'guide_id'}, inplace=True)

    # check if the lengths are the same
    if len(guide_metadata) != len(adata_guide.var):
        print(f"The numbers of sgRNA_ID/guide_id are different: There are {len(adata_guide.var)} in guide anndata, but there are {len(guide_metadata)} in guide metadata.")

    adata_guide.var = adata_guide.var.merge(guide_metadata, left_on='guide_id', right_on ='guide_id', how='inner')
    adata_guide.var[['intended_target_name', 'spacer', 'targeting', 'guide_chr', 'guide_start', 'guide_end', 'intended_target_chr', 'intended_target_start', 'intended_target_end']] = guide_metadata[['intended_target_name', 'spacer', 'targeting', 'guide_chr', 'guide_start', 'guide_end','intended_target_chr', 'intended_target_start', 'intended_target_end']]

    # reset feature_id to var_names (index)
    guide_metadata['feature_id'] = guide_metadata['guide_id'] + "|" + guide_metadata['spacer']
    adata_guide.var_names = guide_metadata['feature_id']

    # adding uns for guide 
    adata_guide.uns['capture_method'] = np.array([capture_method], dtype=object)
    
    # calculate moi
    if moi in ['high', 'low']:
        adata_guide.uns['moi'] = np.array([moi], dtype=object)
    else:
    # Calculate moi if not specified
        binary_matrix = (adata_guide.X > 0).astype(int)
        avg_guides_per_cell = np.mean(np.sum(binary_matrix, axis=1))
        calculated_moi = 'high' if avg_guides_per_cell > 1.5 else 'low'
        adata_guide.uns['moi'] = np.array([calculated_moi], dtype=object)
    
        print(f"Average guides per cell: {avg_guides_per_cell:.2f}")
        print(f"MOI status: {calculated_moi}")

    # adding number of nonzero guides and batch number
    adata_guide.obs['num_expressed_guides'] = (adata_guide.X > 0).sum(axis=1)
    adata_guide.obs['batch_number'] = adata_guide.obs['batch'].factorize()[0] + 1
    adata_guide.obs['total_guide_umis'] = adata_guide.X.sum(axis=1)
    
    ## change in adata_rna; 
    df_gtf['gene_id2'] = df_gtf['gene_id'].str.split('.').str[0]
    df_gtf = df_gtf.drop_duplicates('gene_id2')
    df_gtf_copy = df_gtf.copy()
    df_gtf_copy.set_index('gene_id2', inplace=True)
    # adding gene_start, gene_end, gene_chr
    adata_rna.var = adata_rna.var.join(df_gtf_copy[['seqname', 'start', 'end']].rename(columns={'seqname': 'gene_chr', 
                                                                                            'start': 'gene_start', 
                                                                                            'end': 'gene_end'}))
 
    # rename adata_rna obs
    adata_rna.obs.rename(columns={'n_genes_by_counts': 'n_counts',
                                'pct_counts_mt': 'percent_mito',
                                'n_genes' : 'num_expressed_genes',
                                'total_counts' : 'total_gene_umis'}, inplace=True)
    

    # knee plots
    knee_df = pd.DataFrame({
        'sum': np.array(adata_guide.X.sum(1)).flatten(),
        'barcodes': adata_guide.obs_names.values})
    knee_df = knee_df.sort_values('sum', ascending=False).reset_index(drop=True)
    knee_df['sum_log'] = np.log1p(knee_df['sum'])

    plt.figure(figsize=(8, 5))
    plt.plot(knee_df.index, knee_df['sum_log'], marker='o', linestyle='-', markersize=3)
    plt.xlabel('Barcode Index')
    plt.ylabel('Log of UMI Counts')
    plt.title('Knee Plot')

    # Save knee plot
    if not os.path.exists('figures'):
        os.makedirs('figures')
        print(f"Directory '{'figures'}' created.")
    else:
        print(f"Directory already exists.")

    plt.savefig('figures/knee_plot_guide.png')

    # Find the intersection of barcodes between scRNA and guide data
    intersecting_barcodes = list(set(adata_rna.obs_names)
                                .intersection(adata_guide.obs_names)
                                .intersection(adata_hashing.obs_names))

    mdata = MuData({
        'gene': adata_rna[intersecting_barcodes, :].copy(),
        'guide': adata_guide[intersecting_barcodes, :].copy(),
        'hashing': adata_hashing[intersecting_barcodes, :].copy()
    })

    obs_names = set(mdata.mod['guide'].obs.columns.tolist()) & set(mdata.mod['gene'].obs.columns.tolist()) & set(mdata.mod['hashing'].obs.columns.tolist())
    mdata.obs = mdata.mod['guide'].obs.loc[:, list(obs_names)]

    # Save the MuData object
    mdata.write("mudata.h5mu")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a MuData object from scRNA and guide data.')
    parser.add_argument('adata_rna', type=str, help='Path to the scRNA AnnData file.')
    parser.add_argument('adata_guide', type=str, help='Path to the guide AnnData file.')
    parser.add_argument('adata_hashing', type=str, help ='Path to the hashing AnnData file')
    parser.add_argument('guide_metadata', type=str, help='Path to the guide metadata tsv file.')
    parser.add_argument('gtf', type=str, help='Path to the GTF file.')
    parser.add_argument('moi', default='', help='Multiplicity of infection (MOI) of the screen.')
    parser.add_argument('capture_method', default='', help='Capture Method.')
    
    args = parser.parse_args()
    main(args.adata_rna, args.adata_guide, args.adata_hashing, args.guide_metadata, args.gtf, args.moi, args.capture_method)
