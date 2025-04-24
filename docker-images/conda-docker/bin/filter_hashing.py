#!/usr/bin/env python
import argparse
import anndata as ad

def filter_and_split(hash_file, filtered_rna_file):

    hashing = ad.read_h5ad(hash_file)
    filtered_rna = ad.read_h5ad(filtered_rna_file)

    intersecting_barcodes = list(set(hashing.obs_names).intersection(filtered_rna.obs_names))
    print(f"Number of intersecting barcodes: {len(intersecting_barcodes)}")

    hashing_filtered = hashing[intersecting_barcodes, :]
    batches = hashing_filtered.obs['batch'].unique()
    
    # Save each batch as a separate AnnData file
    for batch in batches:
        adata_subset = hashing_filtered[hashing_filtered.obs['batch'] == batch].copy()
        output_filename = f"{batch}_hashing_filtered.h5ad"
        adata_subset.write(output_filename)
        print(f"Saved {output_filename}")

def main():
    parser = argparse.ArgumentParser(description="Filter and split AnnData by batch")

    parser.add_argument('--hash_file', type=str, required=True, help='Path to the hashing AnnData file')
    parser.add_argument('--rna_file', type=str, required=True, help='Path to the filtered RNA AnnData file')

    args = parser.parse_args()
    filter_and_split(args.hash_file, args.rna_file)

if __name__ == "__main__":
    main()
