#!/usr/bin/env python3
"""
Gene expression mapping QC metrics and visualizations.

This script computes summary statistics and generates distribution plots for
gene expression mapping quality metrics from a MuData (.h5mu) file.

Metrics computed:
  - Total UMI counts per cell (median, mean, std, min, max, q25, q75)
  - Number of expressed genes per cell (median, mean, std, min, max, q25, q75)
  - Mitochondrial percentage per cell (median, mean, std, min, max, q25, q75)
  - Number of cells (total and per batch)

Outputs:
  - Summary metrics TSV (overall and per-batch)
  - Distribution histograms (overall and by batch)
  - Knee plot of UMI counts
"""

import argparse
import logging
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import mudata
import numpy as np
import pandas as pd
from anndata import AnnData

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

sns.set_style("white")


# ---------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------
def compute_column_stats(series: pd.Series) -> dict:
    """Compute summary statistics for a numeric series."""
    return {
        "median": series.median(),
        "mean": series.mean(),
        "std": series.std(),
        "min": series.min(),
        "max": series.max(),
        "q25": series.quantile(0.25),
        "q75": series.quantile(0.75),
    }


def compute_gene_metrics(
    adata: AnnData,
    batch_label: str = "all",
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
) -> dict:
    """
    Compute gene expression mapping metrics for one subset of cells.

    Returns a flat dictionary of metrics suitable for a single row in a DataFrame.
    """
    metrics = {"batch": batch_label, "n_cells": adata.n_obs}

    # UMI counts
    if umi_col in adata.obs.columns:
        stats = compute_column_stats(adata.obs[umi_col])
        for k, v in stats.items():
            metrics[f"umi_{k}"] = v
    else:
        logger.warning(f"Column '{umi_col}' not found in adata.obs")

    # Expressed genes
    if genes_col in adata.obs.columns:
        stats = compute_column_stats(adata.obs[genes_col])
        for k, v in stats.items():
            metrics[f"genes_{k}"] = v
    else:
        logger.warning(f"Column '{genes_col}' not found in adata.obs")

    # Mitochondrial percentage
    if mito_col in adata.obs.columns:
        stats = compute_column_stats(adata.obs[mito_col])
        for k, v in stats.items():
            metrics[f"mito_{k}"] = v
    else:
        logger.warning(f"Column '{mito_col}' not found in adata.obs")

    return metrics


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_knee(
    adata: AnnData,
    outdir: str,
    prefix: str = "gene",
) -> None:
    """Plot knee plot of UMI counts per cell."""
    sums = np.array(adata.X.sum(axis=1)).ravel()
    knee_df = pd.DataFrame({"sum": sums})
    knee_df = knee_df.sort_values("sum", ascending=False).reset_index(drop=True)
    knee_df["sum_log"] = np.log1p(knee_df["sum"])

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(knee_df["sum_log"], knee_df.index, marker="o", linestyle="-", markersize=3)
        ax.set_ylabel("Barcode rank")
        ax.set_xlabel("Log(UMI counts + 1)")
        ax.set_title(f"Knee plot (N={adata.n_obs} cells)")
        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_knee_plot.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_histograms(
    adata: AnnData,
    outdir: str,
    prefix: str = "gene",
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
    bins: int = 100,
) -> None:
    """Plot histograms of gene expression QC metrics (overall, no batch split)."""
    missing = [c for c in [umi_col, genes_col, mito_col] if c not in adata.obs.columns]
    if missing:
        logger.warning(f"Missing columns for histogram plots: {missing}")
        return

    with sns.plotting_context("talk", font_scale=1.2):
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))

        # UMI counts
        sns.histplot(adata.obs[umi_col], bins=bins, ax=axes[0])
        median_val = adata.obs[umi_col].median()
        axes[0].axvline(median_val, color="red", linestyle="--")
        axes[0].text(0.95, 0.95, f"Median: {median_val:.0f}", ha="right", va="top",
                     transform=axes[0].transAxes, fontsize=14)
        axes[0].set_title("UMI counts per cell")
        axes[0].set_xlabel(umi_col)

        # Expressed genes
        sns.histplot(adata.obs[genes_col], bins=bins, ax=axes[1])
        median_val = adata.obs[genes_col].median()
        axes[1].axvline(median_val, color="red", linestyle="--")
        axes[1].text(0.95, 0.95, f"Median: {median_val:.0f}", ha="right", va="top",
                     transform=axes[1].transAxes, fontsize=14)
        axes[1].set_title("Expressed genes per cell")
        axes[1].set_xlabel(genes_col)

        # Mitochondrial %
        sns.histplot(adata.obs[mito_col], bins=bins, ax=axes[2])
        median_val = adata.obs[mito_col].median()
        axes[2].axvline(median_val, color="red", linestyle="--")
        axes[2].text(0.95, 0.95, f"Median: {median_val:.2f}%", ha="right", va="top",
                     transform=axes[2].transAxes, fontsize=14)
        axes[2].set_title("Mitochondrial percentage")
        axes[2].set_xlabel(mito_col)

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_histograms.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_histograms_by_batch(
    adata: AnnData,
    outdir: str,
    prefix: str = "gene",
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
    batch_col: str = "batch",
    bins: int = 100,
) -> None:
    """Plot histograms of gene expression QC metrics split by batch."""
    if batch_col not in adata.obs.columns:
        logger.warning(f"Batch column '{batch_col}' not found; skipping batch-split plots.")
        return

    missing = [c for c in [umi_col, genes_col, mito_col] if c not in adata.obs.columns]
    if missing:
        logger.warning(f"Missing columns for histogram plots: {missing}")
        return

    batches = adata.obs[batch_col].astype("category")
    batch_order = list(batches.cat.categories)
    colors = sns.color_palette("tab10", n_colors=len(batch_order))
    color_map = dict(zip(batch_order, colors))

    def text_y(i: int) -> float:
        return 0.95 - 0.08 * i

    with sns.plotting_context("talk", font_scale=1.0):
        fig, axes = plt.subplots(3, 1, figsize=(10, 10))

        # UMI counts
        sns.histplot(data=adata.obs, x=umi_col, bins=bins, hue=batch_col,
                     palette=color_map, ax=axes[0], element="step",
                     stat="density", common_norm=False)
        for i, batch in enumerate(batch_order):
            median_val = adata.obs.loc[adata.obs[batch_col] == batch, umi_col].median()
            axes[0].axvline(median_val, color=color_map[batch], linestyle="--")
            axes[0].text(0.6, text_y(i), f"{batch}: {median_val:.0f}",
                         ha="right", va="top", transform=axes[0].transAxes,
                         color=color_map[batch], fontsize=12)
        axes[0].set_title("UMI counts per cell")
        axes[0].set_xlabel(umi_col)

        # Expressed genes
        sns.histplot(data=adata.obs, x=genes_col, bins=bins, hue=batch_col,
                     palette=color_map, ax=axes[1], element="step",
                     stat="density", common_norm=False)
        for i, batch in enumerate(batch_order):
            median_val = adata.obs.loc[adata.obs[batch_col] == batch, genes_col].median()
            axes[1].axvline(median_val, color=color_map[batch], linestyle="--")
            axes[1].text(0.8, text_y(i), f"{batch}: {median_val:.0f}",
                         ha="right", va="top", transform=axes[1].transAxes,
                         color=color_map[batch], fontsize=12)
        axes[1].set_title("Expressed genes per cell")
        axes[1].set_xlabel(genes_col)
        if axes[1].legend_ is not None:
            axes[1].legend_.remove()

        # Mitochondrial %
        sns.histplot(data=adata.obs, x=mito_col, bins=bins, hue=batch_col,
                     palette=color_map, ax=axes[2], element="step",
                     stat="density", common_norm=False)
        for i, batch in enumerate(batch_order):
            median_val = adata.obs.loc[adata.obs[batch_col] == batch, mito_col].median()
            axes[2].axvline(median_val, color=color_map[batch], linestyle="--")
            axes[2].text(0.6, text_y(i), f"{batch}: {median_val:.2f}%",
                         ha="right", va="top", transform=axes[2].transAxes,
                         color=color_map[batch], fontsize=12)
        axes[2].set_title("Mitochondrial percentage")
        axes[2].set_xlabel(mito_col)
        if axes[2].legend_ is not None:
            axes[2].legend_.remove()

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_histograms_by_batch.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_cells_per_batch(
    adata: AnnData,
    outdir: str,
    prefix: str = "gene",
    batch_col: str = "batch",
) -> None:
    """Plot bar chart of cells per batch."""
    if batch_col not in adata.obs.columns:
        logger.warning(f"Batch column '{batch_col}' not found; skipping cells-per-batch plot.")
        return

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.countplot(x=batch_col, data=adata.obs, ax=ax, palette="tab10")
        ax.set_title(f"Cells per batch (N={adata.n_obs} total)")
        ax.set_xlabel("Batch")
        ax.set_ylabel("Number of cells")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_cells_per_batch.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def run_gene_mapping_qc(
    input_path: str,
    outdir: str,
    gene_mod_key: str = "gene",
    umi_col: str = "total_gene_umis",
    genes_col: str = "num_expressed_genes",
    mito_col: str = "percent_mito",
    batch_col: str = "batch",
    prefix: str = "gene",
) -> None:
    """
    Run gene expression mapping QC: compute metrics and generate plots.

    Parameters
    ----------
    input_path : str
        Path to MuData (.h5mu) file.
    outdir : str
        Output directory for plots and metrics TSV.
    gene_mod_key : str
        Modality key for gene expression in MuData.
    umi_col, genes_col, mito_col : str
        Column names in adata.obs for QC metrics.
    batch_col : str
        Column name for batch/lane information.
    prefix : str
        Prefix for output filenames.
    """
    os.makedirs(outdir, exist_ok=True)

    # Load data
    logger.info(f"Loading MuData from {input_path}")
    mdata = mudata.read_h5mu(input_path)

    if gene_mod_key not in mdata.mod:
        raise ValueError(f"Modality '{gene_mod_key}' not found in MuData. "
                         f"Available: {list(mdata.mod.keys())}")

    gene = mdata.mod[gene_mod_key]
    logger.info(f"Gene modality: {gene.n_obs} cells, {gene.n_vars} genes")

    # Compute metrics: overall
    metrics_list = []
    metrics_all = compute_gene_metrics(
        gene, batch_label="all",
        umi_col=umi_col, genes_col=genes_col, mito_col=mito_col
    )
    metrics_list.append(metrics_all)

    # Compute metrics: per batch
    if batch_col in gene.obs.columns:
        for batch in sorted(gene.obs[batch_col].unique()):
            adata_batch = gene[gene.obs[batch_col] == batch]
            metrics_batch = compute_gene_metrics(
                adata_batch, batch_label=str(batch),
                umi_col=umi_col, genes_col=genes_col, mito_col=mito_col
            )
            metrics_list.append(metrics_batch)
    else:
        logger.warning(f"Batch column '{batch_col}' not found; skipping per-batch metrics.")

    # Save metrics TSV
    metrics_df = pd.DataFrame(metrics_list)
    metrics_path = os.path.join(outdir, f"{prefix}_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info(f"Saved metrics to {metrics_path}")

    # Generate plots
    plot_knee(gene, outdir, prefix=prefix)
    plot_histograms(gene, outdir, prefix=prefix,
                    umi_col=umi_col, genes_col=genes_col, mito_col=mito_col)
    plot_histograms_by_batch(gene, outdir, prefix=prefix,
                             umi_col=umi_col, genes_col=genes_col,
                             mito_col=mito_col, batch_col=batch_col)
    plot_cells_per_batch(gene, outdir, prefix=prefix, batch_col=batch_col)

    logger.info("Gene mapping QC complete.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Gene expression mapping QC metrics and visualizations."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to MuData (.h5mu) file."
    )
    parser.add_argument(
        "--outdir", "-o", required=True,
        help="Output directory for plots and metrics."
    )
    parser.add_argument(
        "--gene-mod-key", default="gene",
        help="Modality key for gene expression (default: 'gene')."
    )
    parser.add_argument(
        "--umi-col", default="total_gene_umis",
        help="Column name for UMI counts (default: 'total_gene_umis')."
    )
    parser.add_argument(
        "--genes-col", default="num_expressed_genes",
        help="Column name for expressed genes (default: 'num_expressed_genes')."
    )
    parser.add_argument(
        "--mito-col", default="percent_mito",
        help="Column name for mitochondrial % (default: 'percent_mito')."
    )
    parser.add_argument(
        "--batch-col", default="batch",
        help="Column name for batch/lane (default: 'batch')."
    )
    parser.add_argument(
        "--prefix", default="gene",
        help="Prefix for output filenames (default: 'gene')."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_gene_mapping_qc(
        input_path=args.input,
        outdir=args.outdir,
        gene_mod_key=args.gene_mod_key,
        umi_col=args.umi_col,
        genes_col=args.genes_col,
        mito_col=args.mito_col,
        batch_col=args.batch_col,
        prefix=args.prefix,
    )
