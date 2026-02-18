#!/usr/bin/env python3
"""
Chunk MuData for SCEPTRE by splitting only the gene modality.

This script supports auto/off/force chunking modes based on the product
n_cells * n_genes, keeps the full guide modality by default, and filters
pairs_to_test by chunk genes when available.
"""

import argparse
import copy
import os
from math import ceil

import mudata as md
import pandas as pd


def _pairs_to_df(pairs_obj):
    if isinstance(pairs_obj, pd.DataFrame):
        return pairs_obj.copy()
    if isinstance(pairs_obj, dict):
        return pd.DataFrame(pairs_obj)
    return pd.DataFrame(pairs_obj)


def _filter_pairs_by_gene(pairs_df, gene_ids):
    if pairs_df.empty or "gene_id" not in pairs_df.columns:
        return pairs_df
    return pairs_df[pairs_df["gene_id"].isin(set(gene_ids))].copy()


def _copy_global_metadata(source, target):
    target.uns = copy.deepcopy(dict(source.uns))
    if hasattr(source, "obs"):
        target.obs = source.obs.copy()


def _write_single_chunk_passthrough(mdata, mudata_file, output_dir, output_prefix):
    os.makedirs(output_dir, exist_ok=True)
    chunk_path = os.path.join(output_dir, f"{output_prefix}.000.h5mu")

    try:
        src_abs = os.path.abspath(mudata_file)
        if os.path.lexists(chunk_path):
            os.remove(chunk_path)
        os.symlink(src_abs, chunk_path)
    except OSError:
        mdata.write(chunk_path)

    return [chunk_path], [(0, 0, mdata["gene"].n_vars - 1, mdata["gene"].n_vars, "single")]


def chunk_mudata_sceptre(
    mudata_file,
    output_dir,
    chunk_size=4000,
    chunk_mode="auto",
    auto_threshold_entries=2147483647,
    force_chunk=False,
    keep_all_guides=True,
    output_prefix="chunk",
):
    print(f"Loading MuData from {mudata_file}...")
    mdata = md.read(mudata_file)

    if "gene" not in mdata.mod or "guide" not in mdata.mod:
        raise ValueError("MuData must contain both 'gene' and 'guide' modalities")

    n_cells = int(mdata["gene"].n_obs)
    n_genes = int(mdata["gene"].n_vars)
    matrix_entries = int(n_cells) * int(n_genes)

    print(f"Gene modality dimensions: n_cells={n_cells}, n_genes={n_genes}")
    print(f"Gene matrix entries: {matrix_entries}")

    if chunk_mode not in {"auto", "off", "force"}:
        raise ValueError(f"Invalid chunk_mode: {chunk_mode}")

    if force_chunk or chunk_mode == "force":
        should_chunk = True
        mode_used = "force"
    elif chunk_mode == "off":
        should_chunk = False
        mode_used = "off"
    else:
        should_chunk = matrix_entries > int(auto_threshold_entries)
        mode_used = "auto"

    print(
        f"Chunking decision: should_chunk={should_chunk} (mode={mode_used}, threshold={auto_threshold_entries})"
    )

    if not should_chunk or n_genes <= chunk_size:
        files, rows = _write_single_chunk_passthrough(
            mdata=mdata,
            mudata_file=mudata_file,
            output_dir=output_dir,
            output_prefix=output_prefix,
        )
        return files, rows, {
            "mode": mode_used,
            "chunked": False,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "entries": matrix_entries,
            "threshold": int(auto_threshold_entries),
            "chunk_size": int(chunk_size),
            "keep_all_guides": bool(keep_all_guides),
        }

    os.makedirs(output_dir, exist_ok=True)
    n_chunks = ceil(n_genes / int(chunk_size))
    print(f"Splitting {n_genes} genes into {n_chunks} chunks of size {chunk_size}")

    chunk_files = []
    manifest_rows = []

    pairs_df = None
    if "pairs_to_test" in mdata.uns:
        pairs_df = _pairs_to_df(mdata.uns["pairs_to_test"])

    for i in range(n_chunks):
        start = i * int(chunk_size)
        end = min(start + int(chunk_size), n_genes)
        chunk_gene = mdata["gene"][:, start:end].copy()

        if keep_all_guides or pairs_df is None or "guide_id" not in mdata["guide"].var.columns:
            chunk_guide = mdata["guide"].copy()
        else:
            chunk_gene_ids = set(chunk_gene.var_names.tolist())
            filtered_pairs = pairs_df[pairs_df["gene_id"].isin(chunk_gene_ids)]
            keep_guides = set(filtered_pairs.get("guide_id", []))
            chunk_guide = mdata["guide"][:, mdata["guide"].var["guide_id"].isin(keep_guides)].copy()

        modalities = {"gene": chunk_gene, "guide": chunk_guide}

        chunk_mdata = md.MuData(modalities)
        _copy_global_metadata(source=mdata, target=chunk_mdata)

        if pairs_df is not None:
            filtered_pairs = _filter_pairs_by_gene(pairs_df, chunk_gene.var_names.tolist())
            chunk_mdata.uns["pairs_to_test"] = filtered_pairs.to_dict(orient="list")

        chunk_filename = f"{output_prefix}.{i:03d}.h5mu"
        chunk_path = os.path.join(output_dir, chunk_filename)
        chunk_mdata.write(chunk_path)
        chunk_files.append(chunk_path)

        manifest_rows.append((i, start, end - 1, end - start, "chunked"))
        print(f"Saved chunk {i + 1}/{n_chunks}: {chunk_path}")

    summary = {
        "mode": mode_used,
        "chunked": True,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "entries": matrix_entries,
        "threshold": int(auto_threshold_entries),
        "chunk_size": int(chunk_size),
        "keep_all_guides": bool(keep_all_guides),
    }

    return chunk_files, manifest_rows, summary


def write_manifest(manifest_path, manifest_rows, summary, chunk_files):
    rows = []
    for (chunk_id, gene_start, gene_end, n_genes, chunk_mode), chunk_file in zip(manifest_rows, chunk_files):
        rows.append(
            {
                "chunk_id": int(chunk_id),
                "chunk_file": os.path.basename(chunk_file),
                "gene_start": int(gene_start),
                "gene_end": int(gene_end),
                "chunk_gene_count": int(n_genes),
                "chunk_mode": chunk_mode,
                "mode": summary["mode"],
                "chunked": bool(summary["chunked"]),
                "n_cells": int(summary["n_cells"]),
                "n_genes": int(summary["n_genes"]),
                "entries": int(summary["entries"]),
                "threshold": int(summary["threshold"]),
                "chunk_size": int(summary["chunk_size"]),
                "keep_all_guides": bool(summary["keep_all_guides"]),
            }
        )

    manifest_df = pd.DataFrame(rows).sort_values(["chunk_id"]).reset_index(drop=True)
    manifest_df.to_csv(manifest_path, sep="\t", index=False)
    print(f"Wrote manifest: {manifest_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Chunk MuData for SCEPTRE using gene-only chunking",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("mudata_file", help="Input MuData (.h5mu)")
    parser.add_argument("output_dir", help="Output directory for chunks")
    parser.add_argument("--chunk-by", default="gene", choices=["gene"], help="Chunking axis")
    parser.add_argument("--chunk-size", type=int, default=4000, help="Genes per chunk")
    parser.add_argument(
        "--chunk-mode",
        default="auto",
        choices=["auto", "off", "force"],
        help="Chunking mode",
    )
    parser.add_argument(
        "--auto-threshold-entries",
        type=int,
        default=2147483647,
        help="Auto chunk threshold based on n_cells*n_genes",
    )
    parser.add_argument(
        "--force-chunk",
        action="store_true",
        help="Force chunking regardless of mode/threshold",
    )
    parser.add_argument(
        "--keep-all-guides",
        action="store_true",
        help="Keep full guide modality in every chunk",
    )
    parser.add_argument("--output-prefix", default="chunk", help="Chunk file prefix")
    parser.add_argument("--manifest", default="chunk_manifest.tsv", help="Manifest output path")

    args = parser.parse_args()

    if not os.path.exists(args.mudata_file):
        parser.error(f"Input file does not exist: {args.mudata_file}")

    chunk_files, manifest_rows, summary = chunk_mudata_sceptre(
        mudata_file=args.mudata_file,
        output_dir=args.output_dir,
        chunk_size=args.chunk_size,
        chunk_mode=args.chunk_mode,
        auto_threshold_entries=args.auto_threshold_entries,
        force_chunk=args.force_chunk,
        keep_all_guides=args.keep_all_guides,
        output_prefix=args.output_prefix,
    )

    write_manifest(args.manifest, manifest_rows, summary, chunk_files)
    print(f"Created {len(chunk_files)} chunk file(s)")


if __name__ == "__main__":
    raise SystemExit(main())
