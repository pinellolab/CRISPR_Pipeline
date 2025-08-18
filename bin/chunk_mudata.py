#!/usr/bin/env python3
"""
Script to chunk up a mudata file by splitting the gene data into smaller chunks.

This script divides the gene count matrix into chunks of a specified size and
creates separate MuData files for each chunk. It can optionally filter for
specific gene/guide pairs based on mdata.uns["pairs_to_test"].
"""

import argparse
import os
from math import ceil
import pandas as pd
import mudata as md


def chunk_mudata(
    mudata_file,
    output_dir,
    chunk_size=4000,
    test_all_pairs=False,
    output_prefix="chunk",
):
    """
    Split a MuData file into chunks based on the RNA data.

    Parameters:
    -----------
    mudata_file : str
        Path to the input MuData (.h5mu) file
    output_dir : str
        Directory to save the chunked files
    chunk_size : int
        Number of genes per chunk (default: 4000)
    test_all_pairs : bool
        If True, test all pairs; if False, only test pairs in mdata.uns["pairs_to_test"] (default: False)
    output_prefix : str
        Prefix for output chunk files (default: "chunk")

    Returns:
    --------
    list
        List of paths to the created chunk files
    """

    print(f"Loading MuData from {mudata_file}...")
    mdata = md.read(mudata_file)

    # Load pairs_to_test if needed for filtering
    pairs_to_test = None
    if not test_all_pairs:
        if "pairs_to_test" not in mdata.uns:
            raise ValueError(
                "pairs_to_test not found in mdata.uns when test_all_pairs=False"
            )
        print("Loading pairs_to_test from mdata.uns...")
        pairs_to_test = mdata.uns["pairs_to_test"]
        if isinstance(pairs_to_test, dict):
            pairs_to_test = pd.DataFrame(pairs_to_test)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    n_genes = mdata["gene"].shape[1]
    n_chunks = ceil(n_genes / chunk_size)
    print(f"Splitting {n_genes} genes into {n_chunks} chunks of size {chunk_size}")

    chunk_files = []

    for i in range(n_chunks):
        # Define gene range for this chunk
        start = i * chunk_size
        end = min(start + chunk_size, n_genes)
        print(f"Processing chunk {i + 1}/{n_chunks}: genes {start} to {end - 1}")

        # Extract RNA subset for this chunk
        rna_adata_subset = mdata["gene"][:, start:end]

        if not test_all_pairs:
            # Find gRNA elements that target genes in this chunk
            chunk_genes = rna_adata_subset.var_names.tolist()
            chunk_pairs = pairs_to_test[pairs_to_test["gene_id"].isin(chunk_genes)]
            chunk_guides = chunk_pairs["guide_id"].unique()
            print(f"  Found {len(chunk_guides)} unique guides for this chunk")

            # Filter gRNA data to only include relevant guides
            grna_adata_subset = mdata["guide"][
                :, mdata["guide"].var["guide_id"].isin(chunk_guides)
            ]

        else:
            # Use all gRNA data
            grna_adata_subset = mdata["guide"].copy()

        # Create MuData subset
        mdata_subset = md.MuData({"gene": rna_adata_subset, "guide": grna_adata_subset})

        # Copy over pairs_to_test for this chunk if not testing all pairs
        if not test_all_pairs:
            mdata_subset.uns["pairs_to_test"] = chunk_pairs
        else:
            # Copy the original pairs_to_test if it exists
            if "pairs_to_test" in mdata.uns:
                mdata_subset.uns["pairs_to_test"] = mdata.uns["pairs_to_test"]

        print(
            f"  Chunk shape: gene {rna_adata_subset.shape}, guide {grna_adata_subset.shape}"
        )

        # Save chunk
        chunk_filename = f"{output_prefix}.{i:02d}.h5mu"
        chunk_path = os.path.join(output_dir, chunk_filename)
        mdata_subset.write(chunk_path)
        chunk_files.append(chunk_path)
        print(f"  Saved: {chunk_path}")

    print(f"Successfully created {len(chunk_files)} chunk files")
    return chunk_files


def main():
    parser = argparse.ArgumentParser(
        description="Chunk a MuData file by splitting RNA data into smaller chunks",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("mudata_file", help="Path to input MuData (.h5mu) file")

    parser.add_argument("output_dir", help="Directory to save chunked files")

    parser.add_argument(
        "--chunk-size", "-c", type=int, default=4000, help="Number of genes per chunk"
    )

    parser.add_argument(
        "--test-all-pairs",
        action="store_true",
        help="Test all gene/guide pairs instead of only those in pairs_to_test",
    )

    parser.add_argument(
        "--output-prefix", "-p", default="chunk", help="Prefix for output chunk files"
    )

    args = parser.parse_args()

    # Validate arguments
    if not os.path.exists(args.mudata_file):
        parser.error(f"Input file does not exist: {args.mudata_file}")

    # Run chunking
    try:
        chunk_files = chunk_mudata(
            mudata_file=args.mudata_file,
            output_dir=args.output_dir,
            chunk_size=args.chunk_size,
            test_all_pairs=args.test_all_pairs,
            output_prefix=args.output_prefix,
        )

        print("\nChunking completed successfully!")
        print(f"Created {len(chunk_files)} files in {args.output_dir}")

    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
