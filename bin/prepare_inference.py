#!/usr/bin/env python

import argparse
import os
import shutil

import muon as mu
import numpy as np
import pandas as pd

from create_pairs_to_test import generate_pairs_to_test_df
from intended_target_key_utils import annotate_intended_target_groups, enrich_pairs_with_target_metadata


def _close_mudata_handle(mdata):
    file_obj = getattr(mdata, "file", None)
    if file_obj is None:
        return
    close_fn = getattr(file_obj, "close", None)
    if callable(close_fn):
        try:
            close_fn()
        except Exception:
            pass


def _read_pairs_table(path: str) -> pd.DataFrame:
    lower = path.lower()
    if lower.endswith(".parquet"):
        return pd.read_parquet(path)
    if lower.endswith(".tsv") or lower.endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    if lower.endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep=None, engine="python")


def _write_pairs_table(df: pd.DataFrame, prefix: str, output_format: str) -> str:
    if output_format == "parquet":
        parquet_path = f"{prefix}.parquet"
        try:
            df.to_parquet(parquet_path, index=False)
            return parquet_path
        except Exception as err:
            print(
                f"Unable to write parquet pairs table ({err}); falling back to TSV."
            )
    tsv_path = f"{prefix}.tsv"
    df.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


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


def _link_or_copy(src_path: str, dest_path: str):
    if os.path.exists(dest_path):
        os.remove(dest_path)

    src_abs = os.path.abspath(src_path)
    try:
        os.link(src_abs, dest_path)
        return
    except OSError:
        pass

    shutil.copy2(src_abs, dest_path)


def _extract_metadata_backed(mudata_path: str):
    mudata = mu.read_h5mu(mudata_path, backed="r")
    try:
        if "gene" not in mudata.mod or "guide" not in mudata.mod:
            raise KeyError("Mudata file is missing 'gene' or 'guide' modality.")

        guide_var = annotate_intended_target_groups(mudata.mod["guide"].var.copy())
        gene_ids = pd.Index(mudata.mod["gene"].var_names).astype(str)
        return guide_var, gene_ids
    finally:
        _close_mudata_handle(mudata)


def _build_pairs_df(args) -> pd.DataFrame:
    if args.generate_pairs:
        if args.input_gtf is None:
            raise ValueError("--input_gtf is required when --generate_pairs is used.")
        limit = None if args.limit == -1 else args.limit
        return generate_pairs_to_test_df(
            limit=limit,
            mudata_path=args.mudata_path,
            input_gtf=args.input_gtf,
            use_backed=True,
        )

    if args.guide_inference is None:
        raise ValueError(
            "Either --guide_inference or --generate_pairs must be provided."
        )
    return _read_pairs_table(args.guide_inference)


def _subset_mudata_for_cis(mudata_path: str, pairs_to_test: pd.DataFrame) -> mu.MuData:
    mudata = mu.read_h5mu(mudata_path)

    if "gene" not in mudata.mod or "guide" not in mudata.mod:
        raise KeyError("Mudata file is missing 'gene' or 'guide' modality.")

    mudata.mod["guide"].var = annotate_intended_target_groups(mudata.mod["guide"].var)

    tested_guides = pairs_to_test["guide_id"].astype(str).unique()
    tested_genes = pairs_to_test["gene_id"].astype(str).unique()
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

    targeted_cells = np.asarray(grna_subset.X.sum(axis=1)).ravel() > 0
    if not targeted_cells.any():
        raise ValueError("No targeted cells found in the guide modality subset.")

    print(
        f"Keeping {int(targeted_cells.sum())} targeted cells out of {len(targeted_cells)} total cells"
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
    mudata_new.uns["pairs_to_test"] = pairs_to_test.to_dict(orient="list")
    return mudata_new


def main(args):
    print("Preparing guide-gene pairs for inference...")
    pairs_df = _normalize_pairs_schema(_build_pairs_df(args))

    print(f"Loading metadata from {args.mudata_path} using backed reads...")
    guide_var, gene_ids = _extract_metadata_backed(args.mudata_path)

    include_genes = set(pairs_df["gene_name"]).intersection(set(gene_ids))
    include_guides = set(pairs_df["guide_id"]).intersection(set(guide_var["guide_id"].astype(str)))
    print(f"Number of genes in common: {len(include_genes)}")
    print(f"Number of guides in common: {len(include_guides)}")

    subset = pairs_df[
        pairs_df["gene_name"].isin(include_genes)
        & pairs_df["guide_id"].isin(include_guides)
    ].copy()
    if subset.empty:
        raise ValueError(
            "The subset of guide_inference is empty after filtering. Please check your input data."
        )

    subset = enrich_pairs_with_target_metadata(subset, guide_var)
    subset = subset.rename(columns={"gene_name": "gene_id"})
    subset = subset[
        [
            "guide_id",
            "gene_id",
            "intended_target_name",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
            "intended_target_key",
            "pair_type",
        ]
    ].copy()

    pairs_output_path = _write_pairs_table(
        subset,
        prefix=args.output_pairs_prefix,
        output_format=args.pairs_output_format,
    )
    print(f"Wrote sidecar pairs table: {pairs_output_path}")

    if args.subset_for_cis:
        print("Running cis-only mudata subsetting path...")
        mudata_new = _subset_mudata_for_cis(args.mudata_path, subset)
        mudata_new.write(args.output_mudata)
    else:
        print("No cis subsetting requested; passing mudata through without rewrite.")
        _link_or_copy(args.mudata_path, args.output_mudata)

    print(f"Wrote inference mudata artifact: {args.output_mudata}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare pairs_to_test sidecar and mudata for inference."
    )
    parser.add_argument(
        "--mudata_path",
        type=str,
        help="Path to the input MuData file.",
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--guide_inference",
        type=str,
        help="Path to predefined guide-gene pairs table (tsv/csv/parquet).",
    )
    mode.add_argument(
        "--generate_pairs",
        action="store_true",
        help="Generate candidate guide-gene pairs from MuData + GTF.",
    )
    parser.add_argument(
        "--input_gtf",
        type=str,
        default=None,
        help="Input GTF path used only with --generate_pairs.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=1000000,
        help="Distance limit in bp for generated pairs. Use -1 for all genes.",
    )
    parser.add_argument(
        "--subset_for_cis",
        action="store_true",
        help="Subset mudata for cis analysis.",
    )
    parser.add_argument(
        "--pairs_output_format",
        choices=["tsv", "parquet"],
        default="tsv",
        help="Sidecar pairs table format.",
    )
    parser.add_argument(
        "--output_pairs_prefix",
        type=str,
        default="pairs_to_test",
        help="Output prefix for sidecar pairs table.",
    )
    parser.add_argument(
        "--output_mudata",
        type=str,
        default="mudata_inference_input.h5mu",
        help="Output path for prepared MuData file.",
    )
    parser.add_argument(
        "legacy_args",
        nargs="*",
        help=argparse.SUPPRESS,
    )

    args = parser.parse_args()

    # Backward compatibility:
    #   prepare_inference.py <guide_inference> <mudata_path> [--subset_for_cis]
    if args.legacy_args:
        if len(args.legacy_args) == 2 and args.mudata_path is None and args.guide_inference is None and not args.generate_pairs:
            args.guide_inference = args.legacy_args[0]
            args.mudata_path = args.legacy_args[1]
        elif len(args.legacy_args) == 1 and args.generate_pairs and args.mudata_path is None:
            args.mudata_path = args.legacy_args[0]
        else:
            parser.error(
                "Invalid positional arguments. Use either legacy mode "
                "'prepare_inference.py <guide_inference> <mudata_path>' or named arguments."
            )

    if args.mudata_path is None:
        parser.error("--mudata_path is required")
    if args.generate_pairs and args.input_gtf is None:
        parser.error("--input_gtf is required when --generate_pairs is set")
    if not args.generate_pairs and args.guide_inference is None:
        parser.error("Either --guide_inference or --generate_pairs must be provided")

    main(args)
