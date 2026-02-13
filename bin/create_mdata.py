#!/usr/bin/env python

import argparse
import anndata as ad
import pandas as pd
import numpy as np
from muon import MuData
from gtfparse import read_gtf
import matplotlib.pyplot as plt
import os
import sys


def _debug_var_strings(df, label, max_rows=5):
    print(f"[DEBUG] Inspecting {label} var for non-string values")
    cols_to_check = list(df.columns)
    if "spacer" in df.columns:
        cols_to_check = ["spacer"] + [c for c in cols_to_check if c != "spacer"]

    for col in cols_to_check:
        series = df[col]
        dtype = getattr(series.dtype, "name", str(series.dtype))
        if col == "spacer":
            types = series.map(lambda x: type(x).__name__).value_counts(dropna=False)
            print(f"[DEBUG] spacer dtype={dtype} types={types.to_dict()}")

        # Only scan columns that are likely to be written as strings
        is_object = series.dtype == object
        is_categorical = str(series.dtype).startswith("category")
        if not (is_object or is_categorical or col == "spacer"):
            continue

        bad_mask = series.map(lambda x: not (isinstance(x, str) or pd.isna(x)))
        if not bad_mask.any():
            if col == "spacer":
                print(f"[DEBUG] spacer has no non-string values (by isinstance check)")
            continue

        bad_vals = series[bad_mask]
        print(f"[DEBUG] Column '{col}' has {bad_vals.shape[0]} non-string values (dtype={dtype})")
        sample = bad_vals.head(max_rows)
        for idx, val in sample.items():
            print(f"[DEBUG] {label}.var[{idx}]['{col}'] -> {repr(val)} (type={type(val).__name__})")


def _write_debug_files(df, label):
    debug_dir = "debug_create_mdata"
    os.makedirs(debug_dir, exist_ok=True)
    # Write full var table for offline inspection
    var_path = os.path.join(debug_dir, f"{label}_var.tsv")
    df.to_csv(var_path, sep="\t", index=True)

    if "spacer" in df.columns:
        series = df["spacer"]
        all_vals = df.loc[:, ["spacer"]].copy()
        all_vals["spacer_type"] = series.map(lambda x: type(x).__name__)
        all_vals["spacer_isna"] = series.map(lambda x: pd.isna(x))
        all_path = os.path.join(debug_dir, f"{label}_var_spacer_all.tsv")
        all_vals.to_csv(all_path, sep="\t", index=True)

        bad_mask = series.map(lambda x: not (isinstance(x, str) or pd.isna(x)))
        bad_vals = df.loc[bad_mask, ["spacer"]].copy()
        bad_vals["spacer_type"] = series[bad_mask].map(lambda x: type(x).__name__)
        bad_path = os.path.join(debug_dir, f"{label}_var_spacer_nonstring.tsv")
        bad_vals.to_csv(bad_path, sep="\t", index=True)

        type_counts = series.map(lambda x: type(x).__name__).value_counts(dropna=False)
        type_path = os.path.join(debug_dir, f"{label}_var_spacer_type_counts.tsv")
        type_counts.to_csv(type_path, sep="\t", header=["count"])


def _debug_spacer_samples(df, label, max_rows=10):
    if "spacer" not in df.columns:
        print(f"[DEBUG] {label} has no spacer column")
        return
    series = df["spacer"]
    dtype = getattr(series.dtype, "name", str(series.dtype))
    types = series.map(lambda x: type(x).__name__).value_counts(dropna=False)
    na_count = int(pd.isna(series).sum())
    print(f"[DEBUG] {label} spacer dtype={dtype} types={types.to_dict()} na_count={na_count}")
    sample = series.head(max_rows)
    for idx, val in sample.items():
        print(f"[DEBUG] {label}.spacer[{idx}] -> {repr(val)} (type={type(val).__name__})")

def main(adata_rna, adata_guide, guide_metadata, gtf, moi, capture_method, adata_hashing=None, debug_var=False):
    # Load the data
    guide_metadata = pd.read_csv(guide_metadata, sep='\t')
    if debug_var:
        _debug_spacer_samples(guide_metadata, "guide_metadata")
    adata_rna = ad.read_h5ad(adata_rna)
    adata_guide = ad.read_h5ad(adata_guide)
    df_gtf = read_gtf(gtf).to_pandas()

    # Load hashing data if provided
    if 'dummy_hash' in adata_hashing:
        adata_hashing= None
    
    if adata_hashing is not None:
        adata_hashing = ad.read_h5ad(adata_hashing)

    ## change in adata_guide
    # adding var for guide
    adata_guide.var.reset_index(inplace=True)
    # rename gene_id in guide
    adata_guide.var.rename(columns={'gene_id': 'guide_id'}, inplace=True)
    adata_guide.var.rename(columns={'feature_id': 'guide_id'}, inplace=True)


    # check if the lengths are the same
    if len(guide_metadata) != len(adata_guide.var):
        print(f"The numbers of sgRNA_ID/guide_id are different: There are {len(adata_guide.var)} in guide anndata, but there are {len(guide_metadata)} in guide metadata.")

    meta_cols = [
        'intended_target_name', 'spacer', 'targeting', 'guide_chr', 'guide_start',
        'guide_end', 'intended_target_chr', 'intended_target_start', 'intended_target_end'
    ]
    guide_metadata_subset = guide_metadata[['guide_id'] + meta_cols].copy()
    adata_guide.var = adata_guide.var.merge(
        guide_metadata_subset,
        on='guide_id',
        how='left',
        validate='one_to_one'
    )
    if debug_var:
        missing = int(pd.isna(adata_guide.var['spacer']).sum()) if 'spacer' in adata_guide.var.columns else -1
        print(f"[DEBUG] spacer missing after merge: {missing}")
        _debug_spacer_samples(adata_guide.var, "adata_guide.var")

    # reset feature_id to var_names (index)
    assert guide_metadata["guide_id"].is_unique, (
        "guide_id in guide metadata is not unique."
    )
    adata_guide.var_names = guide_metadata["guide_id"]
    
    # adding uns for guide 
    adata_guide.uns['capture_method'] = np.array([capture_method], dtype=object)
    # calculate moi
    if moi in ['high', 'low']:
        adata_guide.uns['moi'] = np.array([moi], dtype=object)
    else:
    # Calculate MOI if not specified
        binary_matrix = (adata_guide.X > 0).astype(int)
        avg_guides_per_cell = np.mean(np.sum(binary_matrix, axis=1))
        calculated_moi = 'high' if avg_guides_per_cell > 1.5 else 'low'
        adata_guide.uns['moi'] = np.array([calculated_moi], dtype=object)
    
        print(f"Average guides per cell: {avg_guides_per_cell:.2f}")
        print(f"MOI status: {calculated_moi}")


    # adding number of nonzero guides and batch number
    adata_guide.obs['num_expressed_guides'] = (adata_guide.X > 0).sum(axis=1)
    adata_guide.obs['total_guide_umis'] = adata_guide.X.sum(axis=1)

    # Add batch_number if hashing data is provided
    if adata_hashing is not None:
        adata_guide.obs['batch_number'] = adata_guide.obs['batch'].factorize()[0] + 1
    
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
                                .intersection(adata_guide.obs_names))

    # Include hashing barcodes if provided
    if adata_hashing is not None:
        intersecting_barcodes = list(set(intersecting_barcodes)
                                   .intersection(adata_hashing.obs_names))

    # Create MuData with conditional hashing modality
    mudata_dict = {
        'gene': adata_rna[intersecting_barcodes, :].copy(),
        'guide': adata_guide[intersecting_barcodes, :].copy()
    }

    if adata_hashing is not None:
        mudata_dict['hashing'] = adata_hashing[intersecting_barcodes, :].copy()

    mdata = MuData(mudata_dict)

    # Get intersection of obs columns across all modalities
    obs_names = set(mdata.mod['guide'].obs.columns.tolist()) & set(mdata.mod['gene'].obs.columns.tolist())
    if adata_hashing is not None:
        obs_names = obs_names & set(mdata.mod['hashing'].obs.columns.tolist())

    mdata.obs = mdata.mod['guide'].obs.loc[:, list(obs_names)]

    # Save the MuData object
    if debug_var:
        _write_debug_files(mdata.mod["guide"].var, "guide")
        _debug_var_strings(mdata.mod["guide"].var, "guide")
        if os.path.isdir("debug_create_mdata"):
            print("[DEBUG] debug_create_mdata contents:", os.listdir("debug_create_mdata"))
        sys.stdout.flush()
        sys.stderr.flush()
    try:
        mdata.write("mudata.h5mu")
    except TypeError:
        if debug_var:
            _debug_var_strings(mdata.mod["guide"].var, "guide")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Create a MuData object from scRNA and guide data.')
    parser.add_argument('adata_rna', type=str, help='Path to the scRNA AnnData file.')
    parser.add_argument('adata_guide', type=str, help='Path to the guide AnnData file.')
    parser.add_argument('guide_metadata', type=str, help='Path to the guide metadata tsv file.')
    parser.add_argument('gtf', type=str, help='Path to the GTF file.')
    parser.add_argument('moi', default='', help='Multiplicity of infection (MOI) of the screen.')
    parser.add_argument('capture_method', default='', help='Capture Method.')
    parser.add_argument('--adata_hashing', type=str, default=None, help='Path to the hashing AnnData file (optional).')
    parser.add_argument('--debug_var', action='store_true', help='Print non-string values in .var when write fails.')

    args = parser.parse_args()
    main(args.adata_rna, args.adata_guide, args.guide_metadata, args.gtf, args.moi, args.capture_method, args.adata_hashing, args.debug_var)
