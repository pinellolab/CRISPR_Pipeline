#!/usr/bin/env python

import argparse
import os
import re

import anndata as ad
import numpy as np
import pandas as pd
from pathlib import Path
from barcode_keys import qualify_barcodes
from count_matrix_utils import describe_matrix, to_sparse_counts


def extract_batch_num(filename):
    match = re.search(r"(.+)_ks_", filename)
    return match.group(1) if match else None


def assert_unique_obs_names(adata, label):
    if adata.obs_names.is_unique:
        return
    duplicates = adata.obs_names[adata.obs_names.duplicated()].unique().tolist()
    raise ValueError(
        f"{label} contains duplicate cell barcodes after concatenation. "
        "The pipeline suffixes barcodes with a shared batch key to keep "
        "measurement sets distinct across modalities. Duplicate examples: "
        + ", ".join(map(str, duplicates[:10]))
    )


def normalize_feature_index(adata, label):
    """Promote an explicit guide_id column to the canonical feature index.

    Older FLASH base-editing mapper outputs used spacer sequences as ``var_names``
    and stored the stable guide identifier in ``adata.var['guide_id']``. AnnData's
    on-disk concatenation does not preserve that auxiliary column by default, so
    normalize those cached outputs before concatenation.
    """
    if "guide_id" not in adata.var.columns:
        return

    guide_ids = adata.var["guide_id"]
    missing = guide_ids.isna() | guide_ids.astype(str).str.strip().eq("")
    if missing.any():
        raise ValueError(
            f"{label} contains {int(missing.sum())} empty guide_id values."
        )

    guide_ids = guide_ids.astype(str)
    duplicates = guide_ids[guide_ids.duplicated(keep=False)].unique().tolist()
    if duplicates:
        raise ValueError(
            f"{label} contains duplicate guide_id values. "
            "Duplicate examples: " + ", ".join(duplicates[:10])
        )

    previous_index = pd.Index(adata.var_names.astype(str))
    if "spacer" not in adata.var.columns and all(
        re.fullmatch(r"[ACGTNacgtn]+", value) for value in previous_index
    ):
        adata.var["spacer"] = previous_index.str.upper().to_numpy()

    adata.var.drop(columns=["guide_id"], inplace=True)
    adata.var_names = pd.Index(guide_ids, name="guide_id")
    print(f"Normalized {label} feature index from guide_id column")


def get_barcode_key(covariates, batch_col, batch_num):
    key_col = "barcode_key" if "barcode_key" in covariates.columns else "concat_batch"
    if key_col not in covariates.columns:
        raise ValueError(
            "parse_covariate.csv is missing both 'barcode_key' and "
            "'concat_batch'. Cannot align cell barcodes across modalities."
        )

    matches = covariates.loc[
        covariates[batch_col].astype(str) == str(batch_num), key_col
    ].dropna().astype(str).unique()
    if len(matches) != 1:
        raise ValueError(
            f"Expected exactly one {key_col} for batch '{batch_num}', "
            f"found {len(matches)}."
        )
    return matches[0]


def apply_batch_suffix(
    adata,
    batch_num,
    barcode_suffix,
    seen_barcodes,
    record_corrected_barcode=False,
):
    obs_names = pd.Index(adata.obs_names.astype(str))
    duplicates = obs_names[obs_names.duplicated()].unique().tolist()
    if duplicates:
        raise ValueError(
            f"Batch {batch_num} contains duplicate cell barcodes. "
            "Duplicate examples: " + ", ".join(map(str, duplicates[:10]))
        )

    suffixed_names = pd.Index(qualify_barcodes(obs_names, barcode_suffix))
    if record_corrected_barcode:
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
            "Batch-qualified cell barcodes are duplicated across mapping outputs. "
            "Duplicate examples: "
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
    parser.add_argument(
        "--normalize-feature-index",
        action="store_true",
        help="Promote an explicit guide_id var column to the canonical feature index.",
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
        if args.normalize_feature_index:
            normalize_feature_index(adata, h5ad_path)
        if all(
            layer in adata.layers
            for layer in ["mature", "nascent", "ambiguous"]
        ):
            # kb-python's nac workflow EM-distributes ambiguous reads across
            # mature/nascent, so these can be genuinely fractional (e.g. 0.5)
            # -- never force them to an integer dtype. Sum in float32 so
            # combining three sparse/dense layers can't silently lose
            # precision to an unexpected int upcast; to_sparse_counts below
            # decides the final on-disk dtype from the actual values, not an
            # assumption about which workflow produced them.
            adata.X = (
                adata.layers["mature"].astype(np.float32)
                + adata.layers["nascent"].astype(np.float32)
                + adata.layers["ambiguous"].astype(np.float32)
            )
            adata.layers["mature"] = to_sparse_counts(adata.layers["mature"])
            adata.layers["nascent"] = to_sparse_counts(adata.layers["nascent"])
            adata.layers["ambiguous"] = to_sparse_counts(adata.layers["ambiguous"])
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

        # Only whole-valued (or already-integer) matrices get narrowed --
        # e.g. mm=False leaves the nac-combined .X genuinely fractional, and
        # it correctly stays float here.
        adata.X = to_sparse_counts(adata.X)
        print(describe_matrix(adata.X, f"{os.path.basename(file_path)} .X"))

        if idx == 0 and adata.var_names.name is not None:
            var_index_name = adata.var_names.name
            print(f"Captured var index name: {var_index_name}")

        batch_num = extract_batch_num(os.path.basename(file_path))
        batch_label = batch_num if batch_num else os.path.basename(file_path)

        barcode_suffix = get_barcode_key(temp, cov_name[0], batch_num)
        apply_batch_suffix(
            adata,
            batch_label,
            barcode_suffix,
            seen_barcodes,
            record_corrected_barcode=bc_replacement,
        )

        if batch_num:
            cov1 = cov_name[0]
            adata.obs[cov1] = batch_num

            adata.obs[cov1] = adata.obs[cov1].astype(str)
            temp[cov1] = temp[cov1].astype(str)

            covariates_for_obs = temp.drop(columns=["barcode_key"], errors="ignore")
            adata.obs = adata.obs.join(covariates_for_obs.set_index(cov1), on=cov1)

        # Save to a permanent file in the temp directory
        processed_file_path = os.path.join(
            args.temp_dir, f"processed_{idx}.h5ad"
        )
        print(f"Saving processed file to {processed_file_path}")
        adata.write_h5ad(processed_file_path)
        processed_files.append(processed_file_path)

    # Use concat_on_disk with the processed file paths
    print(f"Concatenating files: {processed_files}")
    combined_adata = ad.experimental.concat_on_disk(
        [Path(path) for path in processed_files],
        join="outer",
        index_unique=None,
        out_file=Path(args.output)
    )

    final_adata = ad.read_h5ad(args.output)
    if var_index_name:
        print(f"Restoring var index name: {var_index_name}")
        final_adata.var_names.name = var_index_name
    assert_unique_obs_names(final_adata, args.output)
    final_adata.write_h5ad(args.output)

    print(f"Combined AnnData saved to {args.output}")


if __name__ == "__main__":
    main()
