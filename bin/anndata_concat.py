#!/usr/bin/env python

import argparse
import os
import re

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path


def extract_batch_num(filename):
    match = re.search(r"(.+)_ks_", filename)
    return match.group(1) if match else None


def assert_unique_obs_names(adata, label):
    if adata.obs_names.is_unique:
        return
    duplicates = adata.obs_names[adata.obs_names.duplicated()].unique().tolist()
    raise ValueError(
        f"{label} contains duplicate cell barcodes after concatenation. "
        "In barcode replacement mode the pipeline suffixes corrected barcodes "
        "with a shared concat_batch key to avoid modality-specific positional "
        "suffixes. Duplicate examples: "
        + ", ".join(map(str, duplicates[:10]))
    )


def get_concat_batch(covariates, batch_col, batch_num):
    if "concat_batch" not in covariates.columns:
        raise ValueError(
            "parse_covariate.csv is missing 'concat_batch'. "
            "Barcode replacement mode requires a shared cross-modality concat key."
        )

    matches = covariates.loc[
        covariates[batch_col].astype(str) == str(batch_num), "concat_batch"
    ].dropna().astype(str).unique()
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one concat_batch for batch '{batch_num}', "
            f"found {len(matches)}."
        )
    return matches[0]


def apply_replacement_suffix(adata, batch_num, concat_batch, seen_barcodes):
    obs_names = pd.Index(adata.obs_names.astype(str))
    duplicates = obs_names[obs_names.duplicated()].unique().tolist()
    if duplicates:
        raise ValueError(
            f"Batch {batch_num} contains duplicate corrected cell barcodes. "
            "Duplicate examples: " + ", ".join(map(str, duplicates[:10]))
        )

    suffixed_names = pd.Index([f"{barcode}_{concat_batch}" for barcode in obs_names])
    adata.obs["corrected_barcode"] = obs_names
    adata.obs_names = suffixed_names

    repeated = []
    for barcode in suffixed_names:
        previous_batch = seen_barcodes.get(barcode)
        if previous_batch is not None:
            repeated.append((barcode, previous_batch))
        else:
            seen_barcodes[barcode] = batch_num

    if repeated:
        examples = [
            f"{barcode} in {previous_batch} and {batch_num}"
            for barcode, previous_batch in repeated[:10]
        ]
        raise ValueError(
            "Barcode replacement produced corrected cell barcodes that are "
            "duplicated within the same concat batch. Duplicate examples: "
            + "; ".join(examples)
        )


def main():
    parser = argparse.ArgumentParser(
        description="Process AnnData files and add batch information."
    )
    parser.add_argument(
        "file_path_list",
        nargs="+",
        help="File path list for input AnnData files",
    )
    parser.add_argument(
        "parsed_covariate_df",
        type=str,
        help="Path to the parsed covariate DataFrame CSV",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="combined_adata.h5ad",
        help="Output file name (default: combined_adata.h5ad)",
    )
    parser.add_argument(
        "--temp_dir",
        type=str,
        default="./temp_processed",
        help="Directory for processed files",
    )
    parser.add_argument(
        "--bc_replacement",
        type=str,
        help="Barcode replacement specification (true/false)",
        default="false",
    )
    parser.add_argument(
        "--mm",
        type=str,
        help="Multi-mapping specification (true/false)",
        default="false",
    )

    args = parser.parse_args()

    os.makedirs(args.temp_dir, exist_ok=True)

    temp = pd.read_csv(args.parsed_covariate_df)
    cov_name = temp.columns.tolist()
    processed_files = []
    seen_barcodes = {}
    var_index_name = None

    # Sort the file path list by file name
    sorted_file_path_list = sorted(
        args.file_path_list, key=lambda x: os.path.basename(x)
    )
    print(sorted_file_path_list)

    mm = args.mm.lower() in ["true", "1", "yes"]
    bc_replacement = args.bc_replacement.lower() in ["true", "1", "yes"]
    if bc_replacement:
        folder_ = "counts_unfiltered_modified"
    else:
        folder_ = "counts_unfiltered"

    for idx, file_path in enumerate(sorted_file_path_list):
        h5ad_path = os.path.join(file_path, f"{folder_}/adata.h5ad")
        print(f"Processing {h5ad_path}")

        adata = ad.read_h5ad(h5ad_path)
        if all(
            layer in adata.layers
            for layer in ["mature", "nascent", "ambiguous"]
        ):
            adata.X = adata.X.astype(np.float32)
            adata.X = (
                adata.layers["mature"].astype(np.float32)
                + adata.layers["nascent"].astype(np.float32)
                + adata.layers["ambiguous"].astype(np.float32)
            )
            adata.layers["mature"] = adata.layers["mature"].astype(np.float32)
            adata.layers["nascent"] = adata.layers["nascent"].astype(np.float32)
            adata.layers["ambiguous"] = adata.layers["ambiguous"].astype(
                np.float32
            )
            print(
                "Nascent (nac) workflow detected: combining mature, nascent, and ambiguous counts into .X"
            )
        else:
            print("Standard workflow detected: Using existing .X matrix")

        if mm:
            # Round values for downstream compatibility
            if hasattr(adata.X, "toarray"):
                adata.X.data = np.round(adata.X.data)
            else:
                adata.X = np.round(adata.X)

        if idx == 0 and adata.var_names.name is not None:
            var_index_name = adata.var_names.name
            print(f"Captured var index name: {var_index_name}")

        batch_num = extract_batch_num(os.path.basename(file_path))
        batch_label = batch_num if batch_num else os.path.basename(file_path)

        if bc_replacement:
            concat_batch = get_concat_batch(temp, cov_name[0], batch_num)
            apply_replacement_suffix(
                adata, batch_label, concat_batch, seen_barcodes
            )

        if batch_num:
            cov1 = cov_name[0]
            adata.obs[cov1] = batch_num

            adata.obs[cov1] = adata.obs[cov1].astype(str)
            temp[cov1] = temp[cov1].astype(str)

            adata.obs = adata.obs.join(temp.set_index(cov1), on=cov1)

        # Save to a permanent file in the temp directory
        processed_file_path = os.path.join(
            args.temp_dir, f"processed_{idx}.h5ad"
        )
        print(f"Saving processed file to {processed_file_path}")
        adata.write_h5ad(processed_file_path)
        processed_files.append(processed_file_path)

    # Use concat_on_disk with the processed file paths
    print(f"Concatenating files: {processed_files}")
    index_unique = None if bc_replacement else "_"
    combined_adata = ad.experimental.concat_on_disk(
        [Path(path) for path in processed_files],
        join="outer",
        index_unique=index_unique,
        out_file=Path(args.output)
    )

    if var_index_name or bc_replacement:
        final_adata = ad.read_h5ad(args.output)
        if var_index_name:
            print(f"Restoring var index name: {var_index_name}")
            final_adata.var_names.name = var_index_name
        if bc_replacement:
            assert_unique_obs_names(final_adata, args.output)
        final_adata.write_h5ad(args.output)

    print(f"Combined AnnData saved to {args.output}")


if __name__ == "__main__":
    main()
