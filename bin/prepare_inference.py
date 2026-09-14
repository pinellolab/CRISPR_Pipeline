#!/usr/bin/env python

import argparse
import pandas as pd
import mudata as mu
from inference_covariates import materialize_shared_covariates

from analysis_output_formatting import make_h5mu_safe_dataframe
from intended_target_key_utils import (
    annotate_intended_target_groups,
    enrich_pairs_with_target_metadata,
)


def _normalize_pairs_schema(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    if "gene_name" not in out.columns and "gene_id" in out.columns:
        out["gene_name"] = out["gene_id"]

    required_cols = {"guide_id", "gene_name"}
    missing = [col for col in required_cols if col not in out.columns]
    if missing:
        raise ValueError(
            "guide_inference file is missing required columns: " + ", ".join(missing)
        )

    if "pair_type" not in out.columns:
        out["pair_type"] = "discovery"

    out["guide_id"] = out["guide_id"].astype(str)
    out["gene_name"] = out["gene_name"].astype(str)
    return out


def main(guide_inference, mudata_path, subset_for_cis=False):
    print(f"Reading guide inference file from {guide_inference}...")
    guide_inference_df = pd.read_csv(guide_inference)
    guide_inference_df = _normalize_pairs_schema(guide_inference_df)

    print(f"Reading mudata file from {mudata_path}...")
    mudata = mu.read_h5mu(mudata_path)

    if "gene" not in mudata.mod or "guide" not in mudata.mod:
        raise KeyError("Mudata file is missing 'gene' or 'guide' modality.")

    gvar = annotate_intended_target_groups(mudata.mod["guide"].var)
    mudata.mod["guide"].var = gvar

    print(f"Tested pairs dataset contains {len(guide_inference_df)} rows.")
    print(f"Mudata 'gene' modality contains {mudata.mod['gene'].shape[1]} genes.")
    print(f"Mudata 'guide' modality contains {mudata.mod['guide'].shape[1]} guides.")

    gene_var = mudata.mod["gene"].var.index if mudata.mod["gene"].var.index is not None else []
    guide_var = mudata.mod["guide"].var["guide_id"] if mudata.mod["guide"] is not None else []

    include_genes = set(guide_inference_df["gene_name"]).intersection(set(gene_var))
    include_guides = set(guide_inference_df["guide_id"]).intersection(set(guide_var))

    print(f"Number of genes in common: {len(include_genes)}")
    print(f"Number of guides in common: {len(include_guides)}")

    subset = guide_inference_df[
        guide_inference_df["gene_name"].isin(include_genes)
        & guide_inference_df["guide_id"].isin(include_guides)
    ].copy()

    if subset.empty:
        raise ValueError(
            "The subset of guide_inference is empty after filtering. Please check your input data."
        )

    subset = enrich_pairs_with_target_metadata(subset, mudata.mod["guide"].var)
    print(f"Subset contains {len(subset)} rows after filtering and metadata enrichment.")

    key_mapping = {
        "gene_name": "gene_id",
        "guide_id": "guide_id",
        "intended_target_name": "intended_target_name",
        "intended_target_chr": "intended_target_chr",
        "intended_target_start": "intended_target_start",
        "intended_target_end": "intended_target_end",
        "intended_target_key": "intended_target_key",
        "pair_type": "pair_type",
    }
    print("Renaming columns in pairs_to_test...")
    subset = subset.rename(columns=key_mapping)

    # Stored as a DataFrame (not a plain dict of columns) so that low-cardinality
    # string columns (gene/guide ids, target names, chromosomes) can round-trip
    # through .uns as pandas categoricals -- both anndata and R's MuData package
    # preserve this encoding, and it's far more compact on disk than plain strings.
    mudata.uns["pairs_to_test"] = make_h5mu_safe_dataframe(subset)

    if subset_for_cis:
        print("Subsetting mudata for cis analysis...")
        pairs_to_test_df = pd.DataFrame(mudata.uns["pairs_to_test"])

        tested_guides = pairs_to_test_df["guide_id"].astype(str).unique()
        tested_genes = pairs_to_test_df["gene_id"].astype(str).unique()

        print(f"Subsetting to {len(tested_genes)} genes and {len(tested_guides)} guides")

        gene_subset_mask = mudata.mod["gene"].var_names.isin(tested_genes)
        if not gene_subset_mask.any():
            raise ValueError("No tested genes found in gene modality")

        guide_names = mudata.mod["guide"].var["intended_target_name"].astype(str)
        control_mask = guide_names.str.startswith("non-targeting|", na=False) | guide_names.eq("non-targeting")

        guide_subset_mask = mudata.mod["guide"].var["guide_id"].astype(str).isin(tested_guides) | control_mask
        if not guide_subset_mask.any():
            raise ValueError("No tested guides found in guide modality")

        rna_subset = mudata.mod["gene"][:, gene_subset_mask]
        grna_subset = mudata.mod["guide"][:, guide_subset_mask]

        targeted_cells = grna_subset.X.sum(axis=1) > 0
        if not targeted_cells.any():
            raise ValueError("No targeted cells found in the guide modality subset.")

        print(
            f"Keeping {targeted_cells.sum()} targeted cells out of {len(targeted_cells)} total cells"
        )

        mdata_dict = {
            "gene": rna_subset[targeted_cells],
            "guide": grna_subset[targeted_cells],
        }

        for mod in mudata.mod.keys():
            if mod not in mdata_dict:
                mdata_dict[mod] = mudata.mod[mod][targeted_cells]

        mudata_new = mu.MuData(mdata_dict)
        mudata_new.uns = mudata.uns.copy()
        mudata = mudata_new

    output_file = "mudata_inference_input.h5mu"
    print(f"Saving processed mudata to {output_file}...")
    # Both inference methods must condition on the same covariates, so write the
    # agreed set where both look: PerTurbo reads named modality columns, SCEPTRE
    # reads the top-level frame as colData.
    materialize_shared_covariates(mudata)
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
