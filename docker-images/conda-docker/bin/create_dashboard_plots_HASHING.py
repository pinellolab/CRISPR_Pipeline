#!/usr/bin/env python

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import muon as mu
import anndata as ad
import numpy as np 
from umap import UMAP
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy import sparse
import math
import time

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

def plot_cell_counts_HTOs(unfiltered_hashing_demux, save_dir):
    cell_counts = unfiltered_hashing_demux.obs['hto_type_split'].value_counts()
    cell_counts_sort = cell_counts.sort_index()

    plt.figure(figsize=(8, 6))
    cell_counts_sort.plot(kind='bar', color='skyblue')
    plt.yscale('log')
    plt.title('Number of Cells per HTO Type')
    plt.xlabel('HTO Type')
    plt.ylabel('Number of Cells (Log Scale)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plot_path = os.path.join(save_dir, 'cells_per_hto_barplot.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()

def normalize_clr(x):
    if np.sum(x > 0) == 0:
        return x  # Return unchanged if all zeros
    return np.log1p(x / (np.exp(sum(np.log1p(x[x > 0])) / len(x))))

def plot_umap_HTOs(unfiltered_hashing_demux, save_dir, n_pca_components=10):

    print(f"Starting UMAP visualization with HTO data")
    start_time = time.time()
    
    # Set up colors for each unique HTO type
    unique_hto_types = unfiltered_hashing_demux.obs['hto_type_split'].cat.categories
    color_palette = plt.cm.tab20
    colors = [color_palette(i) for i in range(min(len(unique_hto_types), 20))]
    if len(unique_hto_types) > 20:
        colors.extend([color_palette(i % 20) for i in range(20, len(unique_hto_types))])
    hto_color_map = dict(zip(unique_hto_types, colors))
    unfiltered_hashing_demux.obs['hto_color'] = [hto_color_map[hto_type] for hto_type in unfiltered_hashing_demux.obs['hto_type_split']]
    
    # Prepare subplots
    batches = unfiltered_hashing_demux.obs['batch'].unique()
    num_batches = len(batches)
    cols = min(2, num_batches)
    rows = math.ceil(num_batches / cols)
    
    fig, axes = plt.subplots(rows, cols, figsize=(18, 8 * rows))
    if rows * cols > 1:
        axes = axes.flatten()
    else:
        axes = np.array([axes])
    for i, batch in enumerate(batches):
        print(f"Processing {batch} ({i+1}/{num_batches})")
        
        ax = axes[i]
        batch_mask = unfiltered_hashing_demux.obs['batch'] == batch
        batch_cells = np.sum(batch_mask)
    
        if batch_cells < 3:
            print(f"Too few cells in batch {batch}, skipping")
            ax.text(0.5, 0.5, f"Insufficient data: {batch}", ha='center', va='center')
            ax.set_title(f"Batch {batch}: {batch_cells} cells")
            continue
            
        # Get HTO data for this batch
        batch_data = unfiltered_hashing_demux[batch_mask].X
        print(f"Matrix shape: {batch_data.shape[0]} cells and {batch_data.shape[1]} features")
        
        try:
            # Convert sparse matrix to dense if needed
            if sparse.issparse(batch_data):
                data_to_process = batch_data.toarray()
            else:
                data_to_process = batch_data.copy()
            
            # Apply CLR transformation
            transformed_data = np.zeros_like(data_to_process)
            for cell_idx in range(data_to_process.shape[0]):
                transformed_data[cell_idx] = normalize_clr(data_to_process[cell_idx])
            
            # Scale the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(transformed_data)
            
            # Run PCA
            max_components = min(scaled_data.shape[0] - 1, scaled_data.shape[1], 50)
            n_components = min(n_pca_components, max_components)
            
            pca = PCA(n_components=n_components, random_state=42)
            pca_result = pca.fit_transform(scaled_data)
            
            # Run UMAP with lower n_neighbors
            umap_model = UMAP(
                n_neighbors=10,
                min_dist=0.1,
                n_components=2,
                metric='euclidean',
                random_state=42
            )
            
            umap_result = umap_model.fit_transform(pca_result)
            
            # Store UMAP coordinates
            unfiltered_hashing_demux.obs.loc[batch_mask, 'UMAP1'] = umap_result[:, 0]
            unfiltered_hashing_demux.obs.loc[batch_mask, 'UMAP2'] = umap_result[:, 1]
            
            # Create scatter plot
            scatter = ax.scatter(
                unfiltered_hashing_demux.obs.loc[batch_mask, 'UMAP1'],
                unfiltered_hashing_demux.obs.loc[batch_mask, 'UMAP2'],
                c=unfiltered_hashing_demux.obs.loc[batch_mask, 'hto_color'],
                alpha=0.7,
                s=5
            )
            
            # Title and labels
            ax.set_title(f"{batch}, {batch_cells} cells")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            
        except Exception as e:
            print(f"Error processing batch {batch}: {str(e)}")
            ax.text(0.5, 0.5, f"Error: {str(e)}", ha='center', va='center')
            ax.set_title(f"Failed to process: {batch}")
    
    # Remove any empty subplots
    for j in range(num_batches, len(axes)):
        fig.delaxes(axes[j])
    
    # Add a legend
    handles = [plt.Line2D([0], [0], marker='o', color=hto_color_map[k], markersize=8, linestyle='None') 
            for k in unique_hto_types]
    fig.legend(
        handles, 
        [str(k) for k in unique_hto_types], 
        loc='center left', 
        bbox_to_anchor=(1.02, 0.5), 
        title='HTO Type'
    )
    
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    
    # Save the figure
    plot_path = os.path.join(save_dir, 'umap_hto.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_path}")
    plt.close()
    
    print(f"Total time: {time.time() - start_time:.2f} seconds")

def plot_umap_HTOs_singlets(hashing_demux, save_dir, n_pca_components=10):
    
    batch_key = 'batch'
    hto_type_key = 'hto_type_split'
    
    print(f"Starting UMAP visualization with {len(hashing_demux.obs)} cells")
    start_time = time.time()

    # -- 1) Set up colors for each unique HTO type
    unique_hto_types = hashing_demux.obs[hto_type_key].cat.categories
    color_palette = plt.cm.tab20  # Using cm.tab20 is the recommended way
    colors = [color_palette(i) for i in range(min(len(unique_hto_types), 20))]  # limit to 20 colors
    
    # If more than 20 HTO types, cycle through colors
    if len(unique_hto_types) > 20:
        colors.extend([color_palette(i % 20) for i in range(20, len(unique_hto_types))])
        
    hto_color_map = dict(zip(unique_hto_types, colors))
    hashing_demux.obs['hto_color'] = [hto_color_map[hto_type] for hto_type in hashing_demux.obs[hto_type_key]]

    # -- 2) Prepare subplots (rows x columns)
    batches = hashing_demux.obs[batch_key].unique()
    num_batches = len(batches)
    cols = min(2, num_batches) 
    rows = math.ceil(num_batches / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(18, 8 * rows))
    if rows * cols > 1:  # Multiple subplots
        axes = axes.flatten()
    else:  # Single subplot
        axes = np.array([axes])

    # -- 3) For each batch
    for i, batch in enumerate(batches):
        print(f"Processing {batch} ({i+1}/{num_batches})")
        batch_start_time = time.time()
        
        ax = axes[i]
        batch_mask = hashing_demux.obs[batch_key] == batch
        batch_cells = np.sum(batch_mask)
        
        if batch_cells == 0:
            print(f"No cells found for batch {batch}, skipping")
            ax.text(0.5, 0.5, f"No cells in batch: {batch}", ha='center', va='center')
            ax.set_title(f"Empty batch: {batch}")
            continue
            
        print(f"Found {batch_cells} cells in batch {batch}")
        
        # Get the data for this batch - use all features (HTO markers)
        batch_data = hashing_demux[batch_mask].X
        
        # Sanity check for data dimensions
        if batch_data.shape[0] < 3:
            print(f"Too few cells ({batch_data.shape[0]}) in batch {batch} for UMAP, skipping")
            ax.text(0.5, 0.5, f"Too few cells: {batch_data.shape[0]}", ha='center', va='center')
            ax.set_title(f"Insufficient data: {batch}")
            continue
        
        try:
            # Apply CLR transformation first
            print("Applying CLR transformation")
            clr_start_time = time.time()
            
            # Convert sparse matrix to dense if needed
            if sparse.issparse(batch_data):
                print("Converting sparse matrix to dense for CLR")
                data_to_transform = batch_data.toarray()
            else:
                data_to_transform = batch_data
            
            # Apply CLR transformation
            transformed_data = np.zeros_like(data_to_transform)
            for cell_idx in range(data_to_transform.shape[0]):
                transformed_data[cell_idx] = normalize_clr(data_to_transform[cell_idx])
            
            print(f"CLR transformation completed in {time.time() - clr_start_time:.2f} seconds")
            
            # Scale the data after CLR transformation
            print("Scaling CLR-transformed data")
            scaling_start_time = time.time()
            
            # Apply standard scaling 
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(transformed_data)
            print(f"Scaling completed in {time.time() - scaling_start_time:.2f} seconds")
            
            # Run PCA on all features 
            print("Running PCA on all features")
            pca_start_time = time.time()
            
            # Calculate max possible components
            max_components = min(scaled_data.shape[0] - 1, scaled_data.shape[1])
            pca = PCA(n_components=max_components, random_state=42)
            pca_result = pca.fit_transform(scaled_data)
            
            # Select only the top n_pca_components 
            n_dims = min(n_pca_components, pca_result.shape[1])
            pca_subset = pca_result[:, :n_dims]
            print(f"Using top {n_dims} PCA components for UMAP")
            
            # Run UMAP 
            print(f"Running UMAP...")
            umap_start_time = time.time()
            
            umap_model = UMAP(
                n_neighbors=10,
                min_dist=0.1,
                n_components=2,
                metric='euclidean',
                random_state=42
            )
            
            umap_result = umap_model.fit_transform(pca_subset)
            print(f"UMAP completed in {time.time() - umap_start_time:.2f} seconds")
            
            # Create scatter plot
            batch_indices = np.where(batch_mask)[0]
            
            # Map colors to the cells in this batch
            batch_colors = [hashing_demux.obs['hto_color'].iloc[idx] for idx in batch_indices]
            
            # Scatter plot with proper colors
            scatter = ax.scatter(
                umap_result[:, 0],
                umap_result[:, 1],
                c=batch_colors,
                alpha=0.7,
                s=5  
            )
            
            # Title and labels
            ax.set_title(f"{batch}, {batch_cells} cells")
            ax.set_xlabel("UMAP1")
            ax.set_ylabel("UMAP2")
            
        except Exception as e:
            print(f"Error processing batch {batch}: {str(e)}")
            ax.text(0.5, 0.5, f"Error: {str(e)}", ha='center', va='center')
            ax.set_title(f"Failed to process: {batch}")
        
        print(f"Batch {batch} completed in {time.time() - batch_start_time:.2f} seconds")

    # -- 4) Remove any empty subplots if batch count < rows*cols
    for j in range(num_batches, len(axes)):
        fig.delaxes(axes[j])

    # -- 5) Add a legend on the right
    handles = [plt.Line2D([0], [0], marker='o', color=hto_color_map[k], markersize=8, linestyle='None') 
                for k in unique_hto_types]
    fig.legend(
        handles, 
        [str(k) for k in unique_hto_types], 
        loc='center left', 
        bbox_to_anchor=(1.02, 0.5), 
        title='HTO Type'
    )

    plt.tight_layout(rect=[0, 0, 0.85, 1]) 

    # -- 6) Save the figure
    plot_path = os.path.join(save_dir, 'umap_hto_singlets.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {plot_path}")
    plt.close()
    
    print(f"Total time: {time.time() - start_time:.2f} seconds")


def main():
    parser = argparse.ArgumentParser(description="Generate various plots from MuData")
    parser.add_argument('--mudata', required=True, help='Path to the mudata object')
    parser.add_argument('--hashing_demux', required=True, help='Path to the hashing demux anndata file')
    parser.add_argument('--unfiltered_hashing_demux', required=True, help='Path to the unfiltered hashing demux anndata file')
    parser.add_argument('--output_dir', required=True, help='Directory where plots will be saved')

    args = parser.parse_args()

    # Load mudata
    mudata = mu.read_h5mu(args.mudata)
    hashing_demux = ad.read_h5ad(args.hashing_demux)
    unfiltered_hashing_demux = ad.read_h5ad(args.unfiltered_hashing_demux)
    os.makedirs(args.output_dir, exist_ok=True)

    # Generate all plots
    plot_umi_threshold(mudata, args.output_dir)
    plot_detected_genes_threshold(mudata, args.output_dir)
    plot_guides_umi_threshold(mudata, args.output_dir)
    plot_sgRNA_frequencies(mudata, args.output_dir)
    plot_guides_per_cell(mudata, args.output_dir)
    plot_cells_per_guide(mudata, args.output_dir)
    plot_cell_counts_HTOs(unfiltered_hashing_demux, args.output_dir)
    plot_umap_HTOs(unfiltered_hashing_demux, args.output_dir)
    plot_umap_HTOs_singlets(hashing_demux, args.output_dir)

if __name__ == "__main__":
    main()
