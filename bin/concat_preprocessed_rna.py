#!/usr/bin/env python3
"""Concatenate measurement-set-filtered RNA AnnData and apply global gene QC."""

import argparse
import math
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from count_matrix_utils import normalize_sparse_index_dtypes


def validate_feature_index(inputs):
    expected = None
    template_var = None
    for path in inputs:
        current = ad.read_h5ad(path, backed="r")
        index = pd.Index(current.var_names.astype(str))
        if template_var is None:
            template_var = current.var.copy()
        current.file.close()
        if expected is None:
            expected = index
        elif not index.equals(expected):
            raise ValueError(
                "QC-filtered measurement sets do not share the same ordered "
                f"RNA feature index: {path} differs from {inputs[0]}"
            )
    return template_var


def recompute_gene_metrics(adata):
    detected_cells = np.asarray((adata.X > 0).sum(axis=0)).ravel()
    total_counts = np.asarray(adata.X.sum(axis=0)).ravel()
    mean_counts = total_counts / adata.n_obs if adata.n_obs else np.zeros(adata.n_vars)
    adata.var["n_cells_by_counts"] = detected_cells
    adata.var["mean_counts"] = mean_counts
    adata.var["log1p_mean_counts"] = np.log1p(mean_counts)
    adata.var["pct_dropout_by_counts"] = (
        100 * (1 - detected_cells / adata.n_obs) if adata.n_obs else np.nan
    )
    adata.var["total_counts"] = total_counts
    adata.var["log1p_total_counts"] = np.log1p(total_counts)
    return detected_cells


def resolve_min_cells(n_obs, fraction):
    """Resolve the strict fractional support rule used by the pipeline."""
    fraction = float(fraction)
    if not 0 <= fraction < 1:
        raise ValueError("Gene cell-support threshold must be a fraction in [0, 1).")
    return max(1, math.floor(n_obs * fraction) + 1)


def plot_mito_filter(before, after, cutoff, outpath):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    finite = np.asarray(before, dtype=float)
    finite = finite[np.isfinite(finite)]
    bins = np.linspace(0, max(float(np.percentile(finite, 99.5)), cutoff * 1.1, 1), 61)
    for ax, values, title, color in (
        (axes[0], before, f"Before: {len(before):,} cells", "#60a5fa"),
        (axes[1], after, f"After: {len(after):,} cells", "#34d399"),
    ):
        ax.hist(values, bins=bins, color=color, edgecolor="white")
        ax.axvline(cutoff, color="#059669", linestyle=":", label=f"Fixed maximum ({cutoff:g}%)")
        ax.set(title=title, xlabel="Mitochondrial counts (%)", ylabel="Cells")
        ax.legend(frameon=False)
    boxes = axes[2].boxplot(
        [before, after], tick_labels=["Before", "After"], vert=False,
        showfliers=False, patch_artist=True,
    )
    for patch, color in zip(boxes["boxes"], ("#93c5fd", "#6ee7b7")):
        patch.set_facecolor(color)
    axes[2].axvline(cutoff, color="#059669", linestyle=":", label=f"Fixed maximum ({cutoff:g}%)")
    axes[2].set(title="Before/after distribution", xlabel="Mitochondrial counts (%)")
    axes[2].legend(frameon=False)
    fig.suptitle("Post-concatenation mitochondrial filter")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def plot_gene_support(detected_cells, required_cells, n_obs, outpath):
    detected_cells = np.asarray(detected_cells)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].hist(detected_cells, bins=70, color="#60a5fa", edgecolor="white")
    axes[0].axvline(required_cells, color="#7c3aed", linestyle="--", label=f"Required cells = {required_cells:,}")
    axes[0].set(xlabel="Cells detecting gene", ylabel="Genes", title="Absolute gene support")
    fractions = detected_cells / n_obs if n_obs else np.zeros_like(detected_cells, dtype=float)
    axes[1].hist(fractions, bins=70, color="#34d399", edgecolor="white")
    axes[1].axvline(required_cells / n_obs if n_obs else 0, color="#7c3aed", linestyle="--", label="Resolved fraction")
    axes[1].set(xlabel="Fraction of retained cells", ylabel="Genes", title="Fractional gene support")
    for ax in axes:
        ax.legend(frameon=False)
    fig.suptitle("Post-concatenation gene-support filter")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def plot_post_concat_flow(flow, outpath):
    fig, ax = plt.subplots(figsize=(11, 6.5))
    ax.axis("off")
    y = 0.9
    for index, row in flow.iterrows():
        if index:
            ax.annotate("", xy=(0.5, y + 0.03), xytext=(0.5, y + 0.12), arrowprops={"arrowstyle": "->", "color": "#64748b"})
        ax.text(
            0.5, y,
            f"{row['filter']}\n{row['threshold']}\n{int(row['before']):,} → {int(row['after']):,} {row['unit']} ({int(row['removed']):,} removed)",
            ha="center", va="center", fontsize=11,
            bbox={"boxstyle": "round,pad=0.55", "facecolor": "#f8fafc", "edgecolor": "#2563eb"},
        )
        y -= 0.35
    fig.suptitle("Post-concatenation QC flow", fontsize=15)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--output", default="filtered_anndata.h5ad")
    parser.add_argument("--qc-dir", type=Path, default=Path("post_concat_qc"))
    parser.add_argument("--pct-mito", type=float, default=15)
    parser.add_argument("--min-cells-fraction", type=float, default=0.05)
    args = parser.parse_args()
    inputs = sorted(map(Path, args.inputs), key=lambda path: path.name)
    if not inputs:
        parser.error("At least one filtered measurement-set AnnData is required")
    template_var = validate_feature_index(inputs)
    ad.experimental.concat_on_disk(
        inputs,
        join="outer",
        # AnnData 0.11 cannot write pandas Series from merge="first" in
        # concat_on_disk. Feature indices are validated above, so restore the
        # first measurement set's feature metadata explicitly after concat.
        merge=None,
        index_unique=None,
        out_file=Path(args.output),
    )
    combined = ad.read_h5ad(args.output)
    combined.var = template_var.loc[combined.var_names].copy()
    if not combined.obs_names.is_unique:
        raise ValueError("Per-measurement-set QC produced duplicate qualified cell barcodes")
    if "batch" in combined.obs:
        combined.obs["batch_number"] = combined.obs["batch"].factorize()[0] + 1
    combined.X = normalize_sparse_index_dtypes(combined.X)
    if not 0 <= args.pct_mito <= 100:
        raise ValueError("Mitochondrial percentage cutoff must be in [0, 100]")
    if "pct_counts_mt" not in combined.obs:
        raise ValueError("Per-measurement-set QC did not provide pct_counts_mt")
    args.qc_dir.mkdir(parents=True, exist_ok=True)
    cells_before = combined.n_obs
    mito_before = pd.to_numeric(combined.obs["pct_counts_mt"], errors="coerce").to_numpy(dtype=float)
    mito_keep = mito_before < args.pct_mito
    combined = combined[mito_keep].copy()
    mito_after = pd.to_numeric(combined.obs["pct_counts_mt"], errors="coerce").to_numpy(dtype=float)
    plot_mito_filter(
        mito_before, mito_after, args.pct_mito,
        args.qc_dir / "post_concat_mito_before_after.png",
    )

    detected_cells = recompute_gene_metrics(combined)
    required_cells = resolve_min_cells(combined.n_obs, args.min_cells_fraction)
    genes_before = combined.n_vars
    plot_gene_support(
        detected_cells, required_cells, combined.n_obs,
        args.qc_dir / "post_concat_gene_support.png",
    )
    combined = combined[:, detected_cells >= required_cells].copy()
    recompute_gene_metrics(combined)
    combined.write_h5ad(args.output)
    flow = pd.DataFrame([
        {
            "step_order": 1,
            "filter": "QC_pct_mito",
            "threshold": f"pct_counts_mt < {args.pct_mito:g}%",
            "before": cells_before,
            "after": combined.n_obs,
            "removed": cells_before - combined.n_obs,
            "unit": "cells",
        },
        {
            "step_order": 2,
            "filter": "QC_min_cells_per_gene",
            "threshold": f"detected cells >= {required_cells:,} (strictly > {args.min_cells_fraction:g} of retained cells)",
            "before": genes_before,
            "after": combined.n_vars,
            "removed": genes_before - combined.n_vars,
            "unit": "genes",
        },
    ])
    flow.to_csv(args.qc_dir / "post_concat_qc_filter_flow.tsv", sep="\t", index=False)
    plot_post_concat_flow(flow, args.qc_dir / "post_concat_qc_filter_flow.png")
    print(
        f"Concatenated {len(inputs)} QC-filtered measurement sets; retained "
        f"{combined.n_obs}/{cells_before} cells after mitochondrial QC and "
        f"{combined.n_vars}/{genes_before} genes after fractional support QC"
    )


if __name__ == "__main__":
    main()
