#!/usr/bin/env python3
"""
Guide mapping QC metrics and visualizations.

This script computes summary statistics and generates distribution plots for
guide assignment quality metrics from a MuData (.h5mu) file.

Data structures used:
  - guide.obs: Per-cell metadata (one row per cell)
  - guide.var: Per-guide metadata (one row per guide)
  - guide.layers["guide_assignment"]: Binary matrix (cells x guides) indicating assignments

Metrics computed (per batch and overall):

  Per-cell metrics (from guide.obs):
    - total_guide_umis: Total guide UMI counts per cell (median, mean, std, min, max, q25, q75)
    - n_guides_per_cell: Number of guides assigned per cell (mean, std, min, max)
      * Derived from: np.sum(guide.layers["guide_assignment"] > 0, axis=1)

  Assignment summary metrics:
    - n_cells: Total number of cells
    - n_cells_with_guide: Cells with at least 1 guide assigned (n_guides_per_cell > 0)
    - n_cells_exactly_1_guide: Cells with exactly 1 guide assigned (n_guides_per_cell == 1)
    - frac_cells_with_guide: Fraction of cells with any guide assignment
    - max_guides_per_cell: Maximum guides assigned to any single cell

  Per-guide metrics (from guide.var, overall only):
    - n_cells_per_guide: Number of cells each guide is assigned to (median, mean, std)
      * Derived from: np.sum(guide.layers["guide_assignment"] > 0, axis=0)
    - n_guides_total: Total number of guides in the library

Outputs:
  - Summary metrics TSV (overall and per-batch)
  - Distribution histograms (guide UMIs, guides per cell, cells per guide)
  - Knee plot of guide UMI counts
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
from scipy.sparse import issparse

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

sns.set_style("white")


# ---------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------
def compute_guide_assignment_counts(
    guide: AnnData,
    assignment_layer: str = "guide_assignment",
) -> None:
    """
    Compute n_guides_per_cell and n_cells_per_guide from the guide assignment matrix.

    This function adds two columns IN-PLACE:
      - guide.obs["n_guides_per_cell"]: For each cell, count how many guides are assigned
      - guide.var["n_cells_per_guide"]: For each guide, count how many cells it's assigned to

    The guide assignment matrix (guide.layers[assignment_layer]) is expected to be:
      - Shape: (n_cells, n_guides)
      - Values > 0 indicate the guide is assigned to that cell

    Parameters
    ----------
    guide : AnnData
        Guide modality AnnData object.
    assignment_layer : str
        Name of the layer containing the guide assignment matrix.
    """
    if assignment_layer not in guide.layers:
        raise ValueError(
            f"Layer '{assignment_layer}' not found in guide.layers. "
            f"Available layers: {list(guide.layers.keys())}"
        )

    # Get the assignment matrix
    assignment_matrix = guide.layers[assignment_layer]

    # Convert to binary: 1 if assigned (value > 0), 0 otherwise
    # This handles both sparse and dense matrices
    is_assigned = assignment_matrix > 0

    # n_guides_per_cell: sum across columns (axis=1) for each cell
    # Result shape: (n_cells,)
    if issparse(is_assigned):
        n_guides_per_cell = np.asarray(is_assigned.sum(axis=1)).ravel()
    else:
        n_guides_per_cell = np.sum(is_assigned, axis=1)

    # n_cells_per_guide: sum across rows (axis=0) for each guide
    # Result shape: (n_guides,)
    if issparse(is_assigned):
        n_cells_per_guide = np.asarray(is_assigned.sum(axis=0)).ravel()
    else:
        n_cells_per_guide = np.sum(is_assigned, axis=0)

    guide.obs["n_guides_per_cell"] = n_guides_per_cell
    guide.var["n_cells_per_guide"] = n_cells_per_guide

    logger.info(
        f"Computed guide assignment counts from layer '{assignment_layer}': "
        f"n_guides_per_cell (obs), n_cells_per_guide (var)"
    )


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


def compute_guide_metrics(
    guide: AnnData,
    batch_label: str = "all",
    umi_col: str = "total_guide_umis",
    guides_per_cell_col: str = "n_guides_per_cell",
    cells_per_guide_col: str = "n_cells_per_guide",
    include_per_guide_stats: bool = True,
) -> dict:
    """
    Compute guide mapping metrics for one subset of cells.

    Parameters
    ----------
    guide : AnnData
        Guide modality (subset or full).
    batch_label : str
        Label for this batch (e.g., "all" or batch name).
    umi_col : str
        Column in guide.obs with total guide UMIs per cell.
    guides_per_cell_col : str
        Column in guide.obs with number of guides assigned per cell.
    cells_per_guide_col : str
        Column in guide.var with number of cells per guide.
    include_per_guide_stats : bool
        Whether to include per-guide statistics (only meaningful for "all" batch).

    Returns
    -------
    dict
        Flat dictionary of metrics for one row in the output DataFrame.
    """
    metrics = {"batch": batch_label}

    # -----------------------------------------------------------------
    # Cell counts
    # -----------------------------------------------------------------
    n_cells = guide.n_obs
    metrics["n_cells"] = n_cells

    # -----------------------------------------------------------------
    # Guide UMI statistics (from guide.obs[umi_col])
    # -----------------------------------------------------------------
    if umi_col in guide.obs.columns:
        umi_series = guide.obs[umi_col]
        stats = compute_column_stats(umi_series)
        for k, v in stats.items():
            metrics[f"guide_umi_{k}"] = v
    else:
        logger.warning(f"Column '{umi_col}' not found in guide.obs")

    # -----------------------------------------------------------------
    # Guides per cell statistics (from guide.obs[guides_per_cell_col])
    # -----------------------------------------------------------------
    if guides_per_cell_col in guide.obs.columns:
        gpc_series = guide.obs[guides_per_cell_col]

        # Basic statistics
        metrics["guides_per_cell_mean"] = gpc_series.mean()
        metrics["guides_per_cell_std"] = gpc_series.std()
        metrics["guides_per_cell_min"] = gpc_series.min()
        metrics["guides_per_cell_max"] = gpc_series.max()
        metrics["guides_per_cell_median"] = gpc_series.median()

        # Assignment summary counts
        # n_cells_with_guide: count of cells where n_guides_per_cell > 0
        n_cells_with_guide = int((gpc_series > 0).sum())
        metrics["n_cells_with_guide"] = n_cells_with_guide

        # n_cells_exactly_1_guide: count of cells where n_guides_per_cell == 1
        n_cells_exactly_1_guide = int((gpc_series == 1).sum())
        metrics["n_cells_exactly_1_guide"] = n_cells_exactly_1_guide

        # frac_cells_with_guide: fraction of cells with at least one guide
        metrics["frac_cells_with_guide"] = n_cells_with_guide / n_cells if n_cells > 0 else 0.0
    else:
        logger.warning(f"Column '{guides_per_cell_col}' not found in guide.obs")

    # -----------------------------------------------------------------
    # Per-guide statistics (from guide.var[cells_per_guide_col])
    # Only computed for "all" batch since guide.var doesn't change per batch
    # -----------------------------------------------------------------
    if include_per_guide_stats and cells_per_guide_col in guide.var.columns:
        cpg_series = guide.var[cells_per_guide_col]
        metrics["n_guides_total"] = guide.n_vars
        metrics["cells_per_guide_median"] = cpg_series.median()
        metrics["cells_per_guide_mean"] = cpg_series.mean()
        metrics["cells_per_guide_std"] = cpg_series.std()
        metrics["cells_per_guide_min"] = cpg_series.min()
        metrics["cells_per_guide_max"] = cpg_series.max()

    return metrics


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_knee(
    guide: AnnData,
    outdir: str,
    prefix: str = "guide",
    umi_col: str = "total_guide_umis",
) -> None:
    """
    Plot knee plot of guide UMI counts per cell.

    Uses guide.obs[umi_col] rather than summing from X, since guide UMIs
    are typically stored as a pre-computed column.
    """
    if umi_col not in guide.obs.columns:
        logger.warning(f"Column '{umi_col}' not found; skipping knee plot.")
        return

    sums = guide.obs[umi_col].values
    knee_df = pd.DataFrame({"sum": sums})
    knee_df = knee_df.sort_values("sum", ascending=False).reset_index(drop=True)
    knee_df["sum_log"] = np.log1p(knee_df["sum"])

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(10, 5))
        ax.plot(knee_df["sum_log"], knee_df.index, marker="o", linestyle="-", markersize=3)
        ax.set_ylabel("Barcode rank")
        ax.set_xlabel("Log(guide UMI counts + 1)")
        ax.set_title(f"Guide UMI knee plot (N={guide.n_obs} cells)")
        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_knee_plot.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_histograms(
    guide: AnnData,
    outdir: str,
    prefix: str = "guide",
    umi_col: str = "total_guide_umis",
    guides_per_cell_col: str = "n_guides_per_cell",
    cells_per_guide_col: str = "n_cells_per_guide",
    label_col: str = "label",
) -> None:
    """
    Plot histograms of guide QC metrics (overall, no batch split).

    Three panels:
      1. Total guide UMIs per cell (from guide.obs)
      2. Number of guides assigned per cell (from guide.obs)
      3. Number of cells per guide (from guide.var), colored by guide label
    """
    with sns.plotting_context("talk", font_scale=1.2):
        fig, axes = plt.subplots(3, 1, figsize=(12, 10))

        # Panel 1: Guide UMI counts per cell
        if umi_col in guide.obs.columns:
            sns.histplot(guide.obs[umi_col], bins=100, ax=axes[0])
            median_val = guide.obs[umi_col].median()
            axes[0].axvline(median_val, color="red", linestyle="--")
            axes[0].text(0.95, 0.95, f"Median: {median_val:.0f}", ha="right", va="top",
                         transform=axes[0].transAxes, fontsize=14)
            axes[0].set_title("Guide UMI counts per cell")
            axes[0].set_xlabel(umi_col)
        else:
            axes[0].set_title(f"Guide UMI counts (column '{umi_col}' not found)")

        # Panel 2: Guides per cell
        if guides_per_cell_col in guide.obs.columns:
            sns.histplot(guide.obs[guides_per_cell_col], bins=30, ax=axes[1])
            mean_val = guide.obs[guides_per_cell_col].mean()
            axes[1].axvline(mean_val, color="red", linestyle="--")
            axes[1].text(0.95, 0.95, f"Mean: {mean_val:.2f}", ha="right", va="top",
                         transform=axes[1].transAxes, fontsize=14)
            axes[1].set_title("Number of guides assigned per cell")
            axes[1].set_xlabel(guides_per_cell_col)
        else:
            axes[1].set_title(f"Guides per cell (column '{guides_per_cell_col}' not found)")

        # Panel 3: Cells per guide (colored by label if available)
        if cells_per_guide_col in guide.var.columns:
            if label_col in guide.var.columns:
                labels = guide.var[label_col].astype("category")
                label_order = list(labels.cat.categories)
                colors = sns.color_palette("tab10", n_colors=len(label_order))
                color_map = dict(zip(label_order, colors))

                sns.histplot(data=guide.var, x=cells_per_guide_col, bins=30,
                             hue=label_col, palette=color_map, ax=axes[2],
                             element="step", stat="density", common_norm=False)

                # Add median lines for each label
                for i, label in enumerate(label_order):
                    median_val = guide.var.loc[guide.var[label_col] == label, cells_per_guide_col].median()
                    axes[2].axvline(median_val, color=color_map[label], linestyle="--")
            else:
                sns.histplot(guide.var[cells_per_guide_col], bins=30, ax=axes[2])
                median_val = guide.var[cells_per_guide_col].median()
                axes[2].axvline(median_val, color="red", linestyle="--")
                axes[2].text(0.95, 0.95, f"Median: {median_val:.0f}", ha="right", va="top",
                             transform=axes[2].transAxes, fontsize=14)

            axes[2].set_title("Number of cells per guide")
            axes[2].set_xlabel(cells_per_guide_col)
        else:
            axes[2].set_title(f"Cells per guide (column '{cells_per_guide_col}' not found)")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_histograms.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_histograms_by_batch(
    guide: AnnData,
    outdir: str,
    prefix: str = "guide",
    umi_col: str = "total_guide_umis",
    guides_per_cell_col: str = "n_guides_per_cell",
    batch_col: str = "batch",
) -> None:
    """
    Plot histograms of per-cell guide metrics split by batch.

    Two panels:
      1. Total guide UMIs per cell, split by batch
      2. Number of guides assigned per cell, split by batch

    Note: Cells per guide is NOT split by batch because guide.var is shared
    across all cells (a guide's total cell count doesn't make sense per-batch).
    """
    if batch_col not in guide.obs.columns:
        logger.warning(f"Batch column '{batch_col}' not found; skipping batch-split plots.")
        return

    batches = guide.obs[batch_col].astype("category")
    batch_order = list(batches.cat.categories)
    colors = sns.color_palette("tab10", n_colors=len(batch_order))
    color_map = dict(zip(batch_order, colors))

    def text_y(i: int) -> float:
        return 0.95 - 0.08 * i

    with sns.plotting_context("talk", font_scale=1.0):
        fig, axes = plt.subplots(2, 1, figsize=(10, 8))

        # Panel 1: Guide UMIs per cell by batch
        if umi_col in guide.obs.columns:
            sns.histplot(data=guide.obs, x=umi_col, bins=100, hue=batch_col,
                         palette=color_map, ax=axes[0], element="step",
                         stat="density", common_norm=False)
            for i, batch in enumerate(batch_order):
                median_val = guide.obs.loc[guide.obs[batch_col] == batch, umi_col].median()
                axes[0].axvline(median_val, color=color_map[batch], linestyle="--")
                axes[0].text(0.5, text_y(i), f"{batch}: {median_val:.0f}",
                             ha="right", va="top", transform=axes[0].transAxes,
                             color=color_map[batch], fontsize=12)
            axes[0].set_title("Guide UMI counts per cell")
            axes[0].set_xlabel(umi_col)
        else:
            axes[0].set_title(f"Guide UMIs (column '{umi_col}' not found)")

        # Panel 2: Guides per cell by batch (use mean, not median, for count data)
        if guides_per_cell_col in guide.obs.columns:
            sns.histplot(data=guide.obs, x=guides_per_cell_col, bins=30, hue=batch_col,
                         palette=color_map, ax=axes[1], element="step",
                         stat="density", common_norm=False)
            for i, batch in enumerate(batch_order):
                mean_val = guide.obs.loc[guide.obs[batch_col] == batch, guides_per_cell_col].mean()
                axes[1].axvline(mean_val, color=color_map[batch], linestyle="--")
                axes[1].text(0.5, text_y(i), f"{batch}: {mean_val:.2f}",
                             ha="right", va="top", transform=axes[1].transAxes,
                             color=color_map[batch], fontsize=12)
            axes[1].set_title("Number of guides assigned per cell")
            axes[1].set_xlabel(guides_per_cell_col)
            if axes[1].legend_ is not None:
                axes[1].legend_.remove()
        else:
            axes[1].set_title(f"Guides per cell (column '{guides_per_cell_col}' not found)")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_histograms_by_batch.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def run_guide_mapping_qc(
    input_path: str,
    outdir: str,
    guide_mod_key: str = "guide",
    assignment_layer: str = "guide_assignment",
    umi_col: str = "total_guide_umis",
    batch_col: str = "batch",
    label_col: str = "label",
    prefix: str = "guide",
) -> None:
    """
    Run guide mapping QC: compute metrics and generate plots.

    Parameters
    ----------
    input_path : str
        Path to MuData (.h5mu) file.
    outdir : str
        Output directory for plots and metrics TSV.
    guide_mod_key : str
        Modality key for guide data in MuData.
    assignment_layer : str
        Name of layer in guide modality containing guide assignment matrix.
    umi_col : str
        Column name in guide.obs for total guide UMIs per cell.
    batch_col : str
        Column name for batch/lane information.
    label_col : str
        Column name in guide.var for guide labels (e.g., target/NTC).
    prefix : str
        Prefix for output filenames.
    """
    os.makedirs(outdir, exist_ok=True)

    # Load data
    logger.info(f"Loading MuData from {input_path}")
    mdata = mudata.read_h5mu(input_path)

    if guide_mod_key not in mdata.mod:
        raise ValueError(f"Modality '{guide_mod_key}' not found in MuData. "
                         f"Available: {list(mdata.mod.keys())}")

    guide = mdata.mod[guide_mod_key]
    logger.info(f"Guide modality: {guide.n_obs} cells, {guide.n_vars} guides")

    # Compute n_guides_per_cell and n_cells_per_guide from assignment matrix
    compute_guide_assignment_counts(guide, assignment_layer=assignment_layer)

    # Column names for computed metrics
    guides_per_cell_col = "n_guides_per_cell"
    cells_per_guide_col = "n_cells_per_guide"

    # -----------------------------------------------------------------
    # Compute metrics: overall
    # -----------------------------------------------------------------
    metrics_list = []
    metrics_all = compute_guide_metrics(
        guide, batch_label="all",
        umi_col=umi_col,
        guides_per_cell_col=guides_per_cell_col,
        cells_per_guide_col=cells_per_guide_col,
        include_per_guide_stats=True,  # Include guide.var stats for "all"
    )
    metrics_list.append(metrics_all)

    # -----------------------------------------------------------------
    # Compute metrics: per batch
    # -----------------------------------------------------------------
    if batch_col in guide.obs.columns:
        for batch in sorted(guide.obs[batch_col].unique()):
            # Subset to cells in this batch
            guide_batch = guide[guide.obs[batch_col] == batch]

            metrics_batch = compute_guide_metrics(
                guide_batch, batch_label=str(batch),
                umi_col=umi_col,
                guides_per_cell_col=guides_per_cell_col,
                cells_per_guide_col=cells_per_guide_col,
                include_per_guide_stats=False,  # guide.var stats don't make sense per-batch
            )
            metrics_list.append(metrics_batch)
    else:
        logger.warning(f"Batch column '{batch_col}' not found; skipping per-batch metrics.")

    # -----------------------------------------------------------------
    # Save metrics TSV
    # -----------------------------------------------------------------
    metrics_df = pd.DataFrame(metrics_list)
    metrics_path = os.path.join(outdir, f"{prefix}_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info(f"Saved metrics to {metrics_path}")

    # -----------------------------------------------------------------
    # Generate plots
    # -----------------------------------------------------------------
    plot_knee(guide, outdir, prefix=prefix, umi_col=umi_col)
    plot_histograms(
        guide, outdir, prefix=prefix,
        umi_col=umi_col,
        guides_per_cell_col=guides_per_cell_col,
        cells_per_guide_col=cells_per_guide_col,
        label_col=label_col,
    )
    plot_histograms_by_batch(
        guide, outdir, prefix=prefix,
        umi_col=umi_col,
        guides_per_cell_col=guides_per_cell_col,
        batch_col=batch_col,
    )

    logger.info("Guide mapping QC complete.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Guide mapping QC metrics and visualizations."
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
        "--guide-mod-key", default="guide",
        help="Modality key for guide data (default: 'guide')."
    )
    parser.add_argument(
        "--assignment-layer", default="guide_assignment",
        help="Layer name containing guide assignment matrix (default: 'guide_assignment')."
    )
    parser.add_argument(
        "--umi-col", default="total_guide_umis",
        help="Column name for guide UMI counts (default: 'total_guide_umis')."
    )
    parser.add_argument(
        "--batch-col", default="batch",
        help="Column name for batch/lane (default: 'batch')."
    )
    parser.add_argument(
        "--label-col", default="label",
        help="Column name in guide.var for guide labels (default: 'label')."
    )
    parser.add_argument(
        "--prefix", default="guide",
        help="Prefix for output filenames (default: 'guide')."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_guide_mapping_qc(
        input_path=args.input,
        outdir=args.outdir,
        guide_mod_key=args.guide_mod_key,
        assignment_layer=args.assignment_layer,
        umi_col=args.umi_col,
        batch_col=args.batch_col,
        label_col=args.label_col,
        prefix=args.prefix,
    )
