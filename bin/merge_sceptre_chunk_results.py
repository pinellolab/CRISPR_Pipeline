#!/usr/bin/env python3
"""Merge per-chunk SCEPTRE outputs into final tables and MuData."""

import argparse
import os

import mudata as mu
import pandas as pd


def _read_many_tsv(files):
    if not files:
        raise ValueError("No input files supplied")
    dfs = []
    for fp in files:
        if not os.path.exists(fp):
            raise FileNotFoundError(f"Missing file: {fp}")
        dfs.append(pd.read_csv(fp, sep="\t"))
    return pd.concat(dfs, ignore_index=True)


def _assert_unique(df, keys, label, mode="warn"):
    dup_mask = df.duplicated(keys, keep=False)
    if dup_mask.any():
        examples = df.loc[dup_mask, keys].head(10).to_dict(orient="records")
        if mode == "warn":
            print(
                f"Warning: Duplicate {label} keys detected for {keys}; sample duplicates: {examples}"
            )
        elif mode == "error":
            raise ValueError(
                f"Duplicate {label} keys detected for {keys}; sample duplicates: {examples}"
            )


def merge_sceptre_chunk_results(
    per_guide_files,
    per_element_files,
    base_mudata,
    chunk_manifest,
    output_mudata,
    output_per_guide,
    output_per_element,
):
    merged_guide = _read_many_tsv(per_guide_files)
    merged_element = _read_many_tsv(per_element_files)

    required_guide = {"gene_id", "guide_id", "p_value", "log2_fc"}
    required_element = {"gene_id", "intended_target_name", "p_value", "log2_fc"}

    missing_guide = required_guide - set(merged_guide.columns)
    missing_element = required_element - set(merged_element.columns)

    if missing_guide:
        raise ValueError(f"Missing required per-guide columns: {sorted(missing_guide)}")
    if missing_element:
        raise ValueError(f"Missing required per-element columns: {sorted(missing_element)}")

    _assert_unique(merged_guide, ["gene_id", "guide_id"], "per-guide")
    _assert_unique(merged_element, ["gene_id", "intended_target_name"], "per-element")

    merged_guide = merged_guide.sort_values(["gene_id", "guide_id"]).reset_index(drop=True)
    merged_element = (
        merged_element.sort_values(["gene_id", "intended_target_name"]).reset_index(drop=True)
    )

    merged_guide.to_csv(output_per_guide, sep="\t", index=False, compression="gzip")
    merged_element.to_csv(output_per_element, sep="\t", index=False, compression="gzip")

    mdata = mu.read_h5mu(base_mudata)
    mdata.uns["per_guide_results"] = merged_guide
    mdata.uns["per_element_results"] = merged_element

    if chunk_manifest and os.path.exists(chunk_manifest):
        mdata.uns["sceptre_chunk_manifest"] = pd.read_csv(chunk_manifest, sep="\t")

    mdata.write(output_mudata, compression="gzip")


def main():
    parser = argparse.ArgumentParser(
        description="Merge SCEPTRE chunk outputs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--per_guide_files", nargs="+", required=True)
    parser.add_argument("--per_element_files", nargs="+", required=True)
    parser.add_argument("--base_mudata", required=True)
    parser.add_argument("--chunk_manifest", required=True)
    parser.add_argument("--output_mudata", default="inference_mudata.h5mu")
    parser.add_argument("--output_per_guide", default="sceptre_per_guide_output.tsv.gz")
    parser.add_argument("--output_per_element", default="sceptre_per_element_output.tsv.gz")
    args = parser.parse_args()

    merge_sceptre_chunk_results(
        per_guide_files=args.per_guide_files,
        per_element_files=args.per_element_files,
        base_mudata=args.base_mudata,
        chunk_manifest=args.chunk_manifest,
        output_mudata=args.output_mudata,
        output_per_guide=args.output_per_guide,
        output_per_element=args.output_per_element,
    )


if __name__ == "__main__":
    raise SystemExit(main())
