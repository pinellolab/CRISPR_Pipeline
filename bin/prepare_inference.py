#!/usr/bin/env python

import argparse
import pandas as pd
import muon as mu
import numpy as np


def main(guide_inference, mudata_path, subset_for_cis=False):
    # read in files
    print(f"Reading guide inference file from {guide_inference}...")
    guide_inference = pd.read_csv(guide_inference)
    print(f"Reading mudata file from {mudata_path}...")
    mudata = mu.read_h5mu(mudata_path)

    # Check if mudata.mod contains 'gene' and 'guide'
    if "gene" not in mudata.mod or "guide" not in mudata.mod:
        raise KeyError("Mudata file is missing 'gene' or 'guide' modality.")

    # Debugging: Print basic information about the input data
    print(f"Tested pairs dataset contains {len(guide_inference)} rows.")
    print(f"Mudata 'gene' modality contains {mudata.mod['gene'].shape[1]} genes.")
    print(f"Mudata 'guide' modality contains {mudata.mod['guide'].shape[1]} guides.")

    gene_var = (
        mudata.mod["gene"].var.index if mudata.mod["gene"].var.index is not None else []
    )
    guide_var = (
        mudata.mod["guide"].var["guide_id"] if mudata.mod["guide"] is not None else []
    )

    include1 = set(guide_inference["gene_name"]).intersection(set(gene_var))
    include2 = set(guide_inference["guide_id"]).intersection(set(guide_var))

    print(f"Number of genes in common: {len(include1)}")
    print(f"Number of guides in common: {len(include2)}")

    subset = guide_inference[
        guide_inference["gene_name"].isin(include1)
        & guide_inference["guide_id"].isin(include2)
    ]

    # Ensure that subset is not empty
    if subset.empty:
        raise ValueError(
            "The subset of guide_inference is empty after filtering. Please check your input data."
        )

    print(f"Subset contains {len(subset)} rows after filtering.")

    mudata.uns["pairs_to_test"] = subset.to_dict(orient="list")

    # rename keys in uns
    key_mapping = {
        "gene_name": "gene_id",
        "guide_id": "guide_id",
        "intended_target_name": "intended_target_name",
        "pair_type": "pair_type",
    }

    # Set intended_target_name to "non-targeting" for non-targeting guides
    ntc_guide_idx = (
        (~mudata.mod["guide"].var["targeting"])
        | (mudata.mod["guide"].var["type"] == "non-targeting")
        | (mudata.mod["guide"].var["intended_target_name"] == "non-targeting")
    )
    n_ntc_guides = np.sum(ntc_guide_idx)
    print(
        f"{n_ntc_guides} non-targeting guides found. {n_ntc_guides / mudata.mod['guide'].var.shape[0] * 100:.2f}% of total"
    )

    if n_ntc_guides > 0:
        if "non-targeting" not in mudata.mod["guide"].var["intended_target_name"].unique():
            mudata.mod["guide"].var["intended_target_name"].cat.add_categories(["non-targeting"], inplace=True)

        mudata.mod["guide"].var.loc[ntc_guide_idx, "intended_target_name"] = "non-targeting"
    else:
        print("No non-targeting guides found. Make sure you set up the metadata correctly")

    # Ensure pairs_to_test is not None before trying to access items
    if mudata.uns.get("pairs_to_test") is None:
        raise ValueError(
            "'pairs_to_test' in mudata.uns is None, something went wrong when processing the subset."
        )

    # Update the keys in pairs_to_test
    print("Renaming keys in pairs_to_test...")
    mudata.uns["pairs_to_test"] = {
        key_mapping.get(k, k): v for k, v in mudata.uns["pairs_to_test"].items()
    }

    # Subset mudata for cis analysis if requested
    if subset_for_cis:
        print("Subsetting mudata for cis analysis...")
        pairs_to_test_df = pd.DataFrame(mudata.uns["pairs_to_test"])

        # Get unique tested genes and guides
        tested_guides = pairs_to_test_df["guide_id"].unique()
        tested_genes = pairs_to_test_df["gene_id"].unique()

        print(
            f"Subsetting to {len(tested_genes)} genes and {len(tested_guides)} guides"
        )

        # Subset gene modality
        gene_subset_mask = mudata.mod["gene"].var_names.isin(tested_genes)
        if not gene_subset_mask.any():
            raise ValueError("No tested genes found in gene modality")

        # Subset guide modality

        guide_subset_mask = mudata.mod["guide"].var["guide_id"].isin(tested_guides) | (
            mudata.mod["guide"].var["intended_target_name"] == "non-targeting"
        )
        if not guide_subset_mask.any():
            raise ValueError("No tested guides found in guide modality")

        # Create subsets
        rna_subset = mudata.mod["gene"][:, gene_subset_mask]
        grna_subset = mudata.mod["guide"][:, guide_subset_mask]

        # Filter to cells that have at least one targeted guide
        targeted_cells = grna_subset.X.sum(axis=1) > 0
        if not targeted_cells.any():
            raise ValueError("No targeted cells found in the guide modality subset.")

        print(
            f"Keeping {targeted_cells.sum()} targeted cells out of {len(targeted_cells)} total cells"
        )

        # Create new mudata with subsets
        mdata_dict = {
            "gene": rna_subset[targeted_cells],
            "guide": grna_subset[targeted_cells],
        }

        # Copy over any additional modalities
        for mod in mudata.mod.keys():
            if mod not in mdata_dict:
                mdata_dict[mod] = mudata.mod[mod][targeted_cells]

        # Create new mudata object
        mudata_new = mu.MuData(mdata_dict)
        mudata_new.uns = mudata.uns.copy()
        mudata = mudata_new

    # save the mudata
    output_file = "mudata_inference_input.h5mu"
    print(f"Saving processed mudata to {output_file}...")
    mudata.write(output_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process mudata and add guide_inference to uns."
    )
    parser.add_argument(
        "guide_inference", type=str, help="Path to the input inference file"
    )
    parser.add_argument("mudata_path", type=str, help="Path to the input MuData file")
    parser.add_argument(
        "--subset_for_cis",
        action="store_true",
        help="Subset mudata to only genes found in pairs_to_test for cis analysis",
    )

    args = parser.parse_args()
    main(args.guide_inference, args.mudata_path, args.subset_for_cis)
