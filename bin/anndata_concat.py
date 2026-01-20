#!/usr/bin/env python

import anndata as ad
import os
import re
import argparse
import numpy as np
import pandas as pd
import scanpy as sc

from sklearn.mixture import GaussianMixture

def find_threshold_gmm(adata, n_components=2):
    """
    Find a cell threshold based on Gaussian Mixture Model (GMM) fitting to log-transformed UMI counts.
    
    Parameters:
    -----------
    adata: AnnData object
    n_components: int, # Gaussian components (2 for bimodal)
    
    Returns:
    --------
    threshold: UMI count threshold
    """
    
    # Get UMI counts
    if 'total_counts' not in adata.obs.columns:
        sc.pp.calculate_qc_metrics(adata, inplace=True)
    
    umi_counts = adata.obs['total_counts'].values
    log_umis = np.log10(umi_counts + 1).reshape(-1, 1)
    
    # Fit GMM
    gmm = GaussianMixture(n_components=n_components, random_state=42)
    gmm.fit(log_umis)
    
    # Get parameters
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_).flatten()
    weights = gmm.weights_
    
    # Sort by mean (low to high)
    sort_idx = np.argsort(means)
    means = means[sort_idx]
    stds = stds[sort_idx]
    weights = weights[sort_idx]
    
    # Find intersection point between the two Gaussians
    # For two Gaussians, solve: w1*N(μ1,σ1) = w2*N(μ2,σ2)
    from scipy.stats import norm
    
    # Create range between the two means
    x_range = np.linspace(means[0] - 3*stds[0], means[1] + 3*stds[1], 1000)
    
    # Calculate probability densities
    pdf1 = weights[0] * norm.pdf(x_range, means[0], stds[0])
    pdf2 = weights[1] * norm.pdf(x_range, means[1], stds[1])
    
    # Find intersection (where difference is minimal)
    intersection_idx = np.argmin(np.abs(pdf1 - pdf2))
    threshold_log = x_range[intersection_idx]
    threshold = 10**threshold_log - 1
    return threshold

def extract_batch_num(filename):
    match = re.search(r'(.+)_ks_', filename)
    return match.group(1) if match else None

def main():
    parser = argparse.ArgumentParser(description='Process AnnData files and add batch information.')
    parser.add_argument('file_path_list', nargs='+', help='File path list for input AnnData files')
    parser.add_argument('parsed_covariate_df', type=str, help='Path to the parsed covariate DataFrame CSV')
    parser.add_argument('--output', type=str, default='combined_adata.h5ad', help='Output file name (default: combined_adata.h5ad)')
    parser.add_argument('--temp_dir', type=str, default='./temp_processed', help='Directory for processed files')
    parser.add_argument('--modality', type=str, help='Modality')

    args = parser.parse_args()
    
    os.makedirs(args.temp_dir, exist_ok=True)
    
    temp = pd.read_csv(args.parsed_covariate_df)
    cov_name = temp.columns.tolist()
    processed_files = []
    var_index_name = None  

    # Sort the file path list by file name
    sorted_file_path_list = sorted(args.file_path_list, key=lambda x: os.path.basename(x))
    print(sorted_file_path_list)
    
    for idx, file_path in enumerate(sorted_file_path_list):
        h5ad_path = os.path.join(file_path, "counts_unfiltered/adata.h5ad")
        print(f"Processing {h5ad_path}")
        
        adata = ad.read_h5ad(h5ad_path)
        
        if idx == 0 and adata.var_names.name is not None:
            var_index_name = adata.var_names.name
            print(f"Captured var index name: {var_index_name}")
        
        batch_num = extract_batch_num(os.path.basename(file_path))
        
        if batch_num:
            cov1 = cov_name[0]
            adata.obs[cov1] = batch_num
            
            adata.obs[cov1] = adata.obs[cov1].astype(str)
            temp[cov1] = temp[cov1].astype(str)
            
            adata.obs = adata.obs.join(temp.set_index(cov1), on=cov1)

        if args.modality.lower() in ['gex', 'rna']:
            print('Thresholding GEX modality')
            threshold = find_threshold_gmm(adata, n_components=2)
            print(threshold)
            adata = adata[adata.obs['total_counts'] >= threshold, :].copy()
            print(f"Cells retained: {adata.n_obs}")
        
        # Save to a permanent file in the temp directory
        processed_file_path = os.path.join(args.temp_dir, f"processed_{idx}.h5ad")
        print(f"Saving processed file to {processed_file_path}")
        adata.write_h5ad(processed_file_path)
        processed_files.append(processed_file_path)
    
    # Use concat_on_disk with the processed file paths
    print(f"Concatenating files: {processed_files}")
    combined_adata = ad.experimental.concat_on_disk(
        processed_files, 
        join='outer', 
        index_unique="_", 
        out_file=args.output
    )
    
    if var_index_name:
        print(f"Restoring var index name: {var_index_name}")
        final_adata = ad.read_h5ad(args.output)
        final_adata.var_names.name = var_index_name
        final_adata.write_h5ad(args.output)
    
    print(f"Combined AnnData saved to {args.output}")

if __name__ == "__main__":
    main()