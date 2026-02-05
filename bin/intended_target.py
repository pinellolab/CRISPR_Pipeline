#!/usr/bin/env python3
"""
Intended target inference QC metrics and visualizations.

This script evaluates how well guides knock down their intended target genes
using results from trans-regulatory analysis stored in MuData.

NOTE: We use trans results (not cis) because:
  - Trans results include ALL tested guide-gene pairs regardless of genomic distance
  - Cis results only include pairs within a distance threshold
  - Non-targeting guides have no genomic location, so they won't appear in cis results
  - For AUROC/AUPRC evaluation, we need non-targeting guide results as negative controls

Data structures used:
  - mdata.uns["trans_per_guide_results"]: DataFrame with trans test results
      Columns: guide_id, gene_id, sceptre_log2_fc, sceptre_p_value,
               perturbo_log2_fc, perturbo_p_value, ...
  - guide.var: Per-guide metadata
      Key columns: guide_id, intended_target_name, gene_name, label
  - guide.layers["guide_assignment"]: Binary matrix (cells x guides)
  - gene.X or gene.layers[layer]: Gene expression matrix

Processing steps:
  1. Load trans_per_guide_results from mdata.uns
  2. Map gene_id to gene_name using guide.var["intended_target_name"] -> "gene_name"
  3. Extract target name from guide_id (e.g., "CD81#strong" -> "CD81")
  4. Filter to "intended target" tests where gene_name == target_name
  5. Compute knockdown metrics (strong knockdowns, significant tests)
  6. Compute AUROC/AUPRC using targeting vs non-targeting guides
  7. Generate volcano plots, ROC/PR curves, and summary tables

Metrics computed:
  - n_guides_tested: Number of guides with intended target test results
  - n_strong_knockdowns: Guides with log2FC <= threshold (e.g., 60% knockdown)
  - n_significant: Guides with p_value < threshold
  - frac_strong_knockdowns: Fraction of guides with strong knockdown
  - frac_significant: Fraction of significant tests
  - auroc: Area under ROC curve (targeting vs non-targeting discrimination)
  - auprc: Area under precision-recall curve

Outputs:
  - Summary metrics TSV (includes AUROC/AUPRC)
  - Volcano plot (log2FC vs -log10(p_value))
  - ROC and Precision-Recall curves
  - Intended target results TSV (all guides with their knockdown stats)
"""

import argparse
import logging
import os
import sys
from typing import Optional, Literal, Tuple, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import mudata
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from mudata import MuData
from sklearn.metrics import precision_recall_curve, roc_curve, auc

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

sns.set_style("white")


def build_intended_target_map(guide_var: pd.DataFrame) -> pd.Series:
    """
    Build a mapping from intended_target_name (gene ID) to gene_name (symbol).

    Parameters
    ----------
    guide_var : pd.DataFrame
        guide.var DataFrame with columns "intended_target_name" and "gene_name".

    Returns
    -------
    pd.Series
        Index = intended_target_name, values = gene_name.
        Duplicates are dropped (keeps first).
    """
    if "intended_target_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'intended_target_name' column")
    if "gene_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'gene_name' column")

    id_map = guide_var.set_index("intended_target_name")["gene_name"]
    id_map = id_map[~id_map.index.duplicated(keep="first")]
    return id_map


# ---------------------------------------------------------------------
# Core QC functions
# ---------------------------------------------------------------------
def load_inference_results(
    mdata: MuData,
    results_key: str = "trans_per_guide_results",
) -> pd.DataFrame:
    """
    Load perturbation inference results from MuData.uns.

    Parameters
    ----------
    mdata : MuData
        MuData object.
    results_key : str
        Key in mdata.uns containing results DataFrame.
        Default is "trans_per_guide_results" which includes all guide-gene pairs.

    Returns
    -------
    pd.DataFrame
        Test results with columns like guide_id, gene_id, log2_fc, p_value.
    """
    if results_key not in mdata.uns:
        raise ValueError(
            f"Results key '{results_key}' not found in mdata.uns. "
            f"Available keys: {list(mdata.uns.keys())}"
        )

    results = mdata.uns[results_key]
    if isinstance(results, pd.DataFrame):
        return results.copy()
    if isinstance(results, dict):
        return pd.DataFrame(results)
    try:
        return pd.DataFrame(results)
    except Exception as exc:
        raise ValueError(
            f"Results at key '{results_key}' could not be converted to a DataFrame."
        ) from exc


def resolve_results_key(mdata: MuData, results_key: str) -> str:
    """Resolve an inference results key, with auto-detection fallback."""
    available = set(mdata.uns.keys())
    if results_key != "auto" and results_key in available:
        return results_key

    if results_key != "auto" and results_key not in available:
        logger.warning(
            f"Requested results key '{results_key}' not found. Falling back to auto-detection."
        )

    # Prefer trans results, then generic per-guide results
    candidates = [
        "trans_per_guide_results",
        "per_guide_results",
        "cis_per_guide_results",
        "trans_test_results",
        "test_results",
    ]
    for key in candidates:
        if key in available:
            logger.info(f"Using results key '{key}' for intended target QC")
            return key

    raise ValueError(
        f"No compatible inference results found in mdata.uns. Available keys: {list(available)}"
    )


def resolve_metric_columns(
    results: pd.DataFrame,
    log2fc_col: str,
    pvalue_col: str,
) -> Tuple[str, str]:
    """Resolve log2FC and p-value column names, with auto-detection fallback."""
    log2fc_candidates = ["log2_fc", "perturbo_log2_fc", "sceptre_log2_fc"]
    pvalue_candidates = ["p_value", "perturbo_p_value", "sceptre_p_value", "q_value"]

    def pick_col(label: str, requested: str, candidates: List[str]) -> str:
        if requested != "auto" and requested in results.columns:
            return requested
        if requested != "auto":
            logger.warning(
                f"Requested {label} column '{requested}' not found. Falling back to auto-detection."
            )
        for cand in candidates:
            if cand in results.columns:
                logger.info(f"Using '{cand}' for {label}")
                return cand
        raise ValueError(
            f"Could not resolve {label} column. Tried: {candidates}. "
            f"Available columns: {results.columns.tolist()}"
        )

    return pick_col("log2fc", log2fc_col, log2fc_candidates), pick_col("pvalue", pvalue_col, pvalue_candidates)


def filter_to_intended_targets(
    trans_results: pd.DataFrame,
    guide_var: pd.DataFrame,
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
) -> pd.DataFrame:
    """
    Filter trans results to only intended target tests (guide targets its intended gene).

    Processing steps:
      1. Map gene_id to gene_name using guide.var
      2. Extract target_name from guide_id (e.g., "CD81#strong" -> "CD81")
      3. Keep rows where gene_name == target_name

    Parameters
    ----------
    trans_results : pd.DataFrame
        Trans results with guide_id, gene_id columns.
    guide_var : pd.DataFrame
        guide.var with intended_target_name and gene_name columns.
    log2fc_col : str
        Column name for log2 fold change.
    pvalue_col : str
        Column name for p-value.

    Returns
    -------
    pd.DataFrame
        Filtered results with added columns: gene_name, target_name.
    """
    # Build gene ID -> gene name mapping
    id_map = build_intended_target_map(guide_var)

    # Add gene_name column by mapping gene_id
    results = trans_results.copy()
    results["gene_name"] = results["gene_id"].map(id_map)

    # Extract target_name from guide_id (format: "TARGET#variant" -> "TARGET")
    results["target_name"] = results["guide_id"].str.split("#").str[0]

    # Filter to intended targets only
    intended = results[results["gene_name"] == results["target_name"]].copy()

    # Deduplicate by guide_id + gene_id
    intended = intended.drop_duplicates(subset=["guide_id", "gene_id"])

    # Merge with full guide metadata
    guide_meta = guide_var.reset_index(drop=True)
    intended = guide_meta.merge(
        intended[["guide_id", "gene_id", log2fc_col, pvalue_col]],
        on="guide_id",
        how="left",
    )

    return intended


def compute_knockdown_metrics(
    intended_results: pd.DataFrame,
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
    fc_threshold: float = 0.4,
    pval_threshold: float = 0.05,
) -> dict:
    """
    Compute summary metrics for intended target knockdown efficiency.

    Parameters
    ----------
    intended_results : pd.DataFrame
        Filtered cis results for intended targets.
    log2fc_col : str
        Column with log2 fold change values.
    pvalue_col : str
        Column with p-values.
    fc_threshold : float
        Fold change threshold for "strong" knockdown.
        Default 0.4 means 60% knockdown (log2(0.4) = -1.32).
    pval_threshold : float
        P-value threshold for significance.

    Returns
    -------
    dict
        Dictionary of metrics.
    """
    # Filter to rows with valid results
    valid = intended_results.dropna(subset=[log2fc_col, pvalue_col])

    n_guides_tested = len(valid)
    log2fc_threshold = np.log2(fc_threshold)

    # Strong knockdowns: log2FC <= threshold (negative = knockdown)
    strong_knockdowns = valid[valid[log2fc_col] <= log2fc_threshold]
    n_strong = len(strong_knockdowns)

    # Significant tests
    significant = valid[valid[pvalue_col] < pval_threshold]
    n_sig = len(significant)

    # Both strong AND significant
    strong_and_sig = valid[
        (valid[log2fc_col] <= log2fc_threshold) & (valid[pvalue_col] < pval_threshold)
    ]
    n_strong_and_sig = len(strong_and_sig)

    metrics = {
        "n_guides_total": len(intended_results),
        "n_guides_tested": n_guides_tested,
        "fc_threshold": fc_threshold,
        "log2fc_threshold": log2fc_threshold,
        "pval_threshold": pval_threshold,
        "n_strong_knockdowns": n_strong,
        "n_significant": n_sig,
        "n_strong_and_significant": n_strong_and_sig,
        "frac_strong_knockdowns": n_strong / n_guides_tested if n_guides_tested > 0 else 0,
        "frac_significant": n_sig / n_guides_tested if n_guides_tested > 0 else 0,
        "frac_strong_and_significant": n_strong_and_sig / n_guides_tested if n_guides_tested > 0 else 0,
        "median_log2fc": valid[log2fc_col].median() if n_guides_tested > 0 else np.nan,
        "mean_log2fc": valid[log2fc_col].mean() if n_guides_tested > 0 else np.nan,
    }

    return metrics


# ---------------------------------------------------------------------
# AUROC/AUPRC Evaluation
# ---------------------------------------------------------------------
def build_evaluation_table(
    trans_results: pd.DataFrame,
    guide_var: pd.DataFrame,
    pvalue_col: str = "p_value",
    non_targeting_label: str = "non_targeting",
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Build evaluation table for AUROC/AUPRC calculation.

    This follows the approach from CRISPR_Pipeline's evaluate_controls.py:
      - Positive controls (direct_target=1): Tests where guide targets its intended gene
      - Negative controls (direct_target=0): Non-targeting guides tested against intended
        target genes (provides matched controls with no expected effect)

    The dataset is balanced by downsampling negative controls to match the number
    of positive controls.

    NOTE: We use trans results (not cis) because cis results only include guide-gene
    pairs within a genomic distance threshold. Non-targeting guides have no genomic
    location, so they won't appear in cis results. Trans results include all tested
    guide-gene pairs regardless of genomic distance.

    Parameters
    ----------
    trans_results : pd.DataFrame
        Trans test results containing all guide-gene pairs tested.
        Must have columns: guide_id, gene_id, and the p-value column.
    guide_var : pd.DataFrame
        guide.var with columns: guide_id, intended_target_name, gene_name, label.
    pvalue_col : str
        Column containing p-values.
    non_targeting_label : str
        Value in guide.var["label"] indicating non-targeting guides.
    random_state : int
        Random seed for balanced sampling.

    Returns
    -------
    pd.DataFrame
        Evaluation table with columns: guide_id, gene_id, p_value, direct_target.
    """
    # Build gene ID -> gene name mapping
    id_map = build_intended_target_map(guide_var)

    # Prepare results with gene_name and target_name
    results = trans_results.copy()
    results["gene_name"] = results["gene_id"].map(id_map)
    results["target_name"] = results["guide_id"].str.split("#").str[0]

    # Merge guide metadata to get label
    # Handle case where guide_id is both index and column
    guide_meta = guide_var.reset_index(drop=True)[["guide_id", "label"]].drop_duplicates()
    results = results.merge(guide_meta, on="guide_id", how="left")

    # Get set of all intended target genes (gene symbols)
    intended_target_genes = set(guide_var["gene_name"].dropna().unique())

    # -------------------------------------------------------------------------
    # Positive controls: guide tests its own intended target gene
    # gene_name == target_name AND label != non_targeting (targeting guides)
    # -------------------------------------------------------------------------
    positives = results[
        (results["gene_name"] == results["target_name"]) &
        (results["label"] != non_targeting_label) &
        (results[pvalue_col].notna())
    ].copy()
    positives["direct_target"] = 1

    # -------------------------------------------------------------------------
    # Negative controls: non-targeting guides tested against intended target genes
    # label == non_targeting AND gene_name is in intended_target_genes
    # -------------------------------------------------------------------------
    negatives = results[
        (results["label"] == non_targeting_label) &
        (results["gene_name"].isin(intended_target_genes)) &
        (results[pvalue_col].notna())
    ].copy()
    negatives["direct_target"] = 0

    # Balance classes by downsampling negatives
    n_pos = len(positives)
    n_neg = len(negatives)

    if n_neg > n_pos and n_pos > 0:
        negatives = negatives.sample(n=n_pos, random_state=random_state)
    elif n_neg == 0 or n_pos == 0:
        logger.warning(f"Cannot build evaluation table: {n_pos} positives, {n_neg} negatives")
        return pd.DataFrame()

    # Combine
    eval_table = pd.concat([positives, negatives], ignore_index=True)
    eval_table = eval_table[["guide_id", "gene_id", pvalue_col, "direct_target"]].copy()
    eval_table = eval_table.rename(columns={pvalue_col: "p_value"})

    logger.info(f"Built evaluation table: {n_pos} positives, {len(negatives)} negatives")

    return eval_table


def compute_auroc_auprc(
    eval_table: pd.DataFrame,
) -> Tuple[float, float, dict]:
    """
    Compute AUROC and AUPRC from evaluation table.

    Uses sklearn's roc_curve and precision_recall_curve with 1-p_value as the
    prediction score (lower p-values indicate positive class).

    Parameters
    ----------
    eval_table : pd.DataFrame
        Must have columns: p_value, direct_target.

    Returns
    -------
    auroc : float
        Area under ROC curve.
    auprc : float
        Area under precision-recall curve.
    curves : dict
        Dictionary containing FPR, TPR, precision, recall arrays for plotting.
    """
    if len(eval_table) == 0:
        return np.nan, np.nan, {}

    true_labels = eval_table["direct_target"].values
    pred_scores = 1 - eval_table["p_value"].values  # Invert: low p-value = high score

    # ROC curve and AUROC
    fpr, tpr, _ = roc_curve(true_labels, pred_scores)
    auroc = auc(fpr, tpr)

    # Precision-recall curve and AUPRC
    precision, recall, _ = precision_recall_curve(true_labels, pred_scores)
    auprc = auc(recall, precision)

    curves = {
        "fpr": fpr,
        "tpr": tpr,
        "precision": precision,
        "recall": recall,
    }

    return auroc, auprc, curves


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------
def plot_roc_pr_curves(
    curves: dict,
    auroc: float,
    auprc: float,
    outdir: str,
    prefix: str = "intended_target",
) -> None:
    """
    Plot ROC and Precision-Recall curves side by side.

    Parameters
    ----------
    curves : dict
        Dictionary with fpr, tpr, precision, recall arrays from compute_auroc_auprc.
    auroc : float
        AUROC value for annotation.
    auprc : float
        AUPRC value for annotation.
    outdir : str
        Output directory.
    prefix : str
        Filename prefix.
    """
    if not curves:
        logger.warning("No curve data available for ROC/PR plot")
        return

    with sns.plotting_context("talk"):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # ROC curve
        ax = axes[0]
        ax.plot(curves["fpr"], curves["tpr"], color="#3498db", linewidth=2)
        ax.plot([0, 1], [0, 1], "k--", alpha=0.5, linewidth=1)
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_title(f"ROC Curve (AUROC = {auroc:.3f})")
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1.05])
        ax.set_aspect("equal")

        # PR curve
        ax = axes[1]
        ax.plot(curves["recall"], curves["precision"], color="#e74c3c", linewidth=2)
        ax.axhline(0.5, color="k", linestyle="--", alpha=0.5, linewidth=1)
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_title(f"Precision-Recall Curve (AUPRC = {auprc:.3f})")
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1.05])
        ax.set_aspect("equal")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_roc_pr_curves.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_volcano(
    intended_results: pd.DataFrame,
    outdir: str,
    prefix: str = "intended_target",
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
    label_col: str = "label",
    pval_threshold: float = 0.05,
    n_label: int = 10,
) -> None:
    """
    Plot volcano plot of intended target knockdown results.

    X-axis: log2 fold change
    Y-axis: -log10(p-value)
    Points colored by guide label (e.g., positive_control, target, non_targeting)
    Top N most significant guides are labeled.

    Parameters
    ----------
    intended_results : pd.DataFrame
        Filtered results for intended targets.
    outdir : str
        Output directory.
    prefix : str
        Prefix for output filename.
    log2fc_col, pvalue_col : str
        Column names for fold change and p-value.
    label_col : str
        Column for coloring points.
    pval_threshold : float
        Significance threshold for horizontal line.
    n_label : int
        Number of top guides to label.
    """
    plot_data = intended_results.dropna(subset=[log2fc_col, pvalue_col]).copy()
    if len(plot_data) == 0:
        logger.warning("No valid data for volcano plot")
        return

    # Compute -log10(p-value), adding small epsilon to avoid log(0)
    plot_data["neg_log10_pval"] = -np.log10(plot_data[pvalue_col] + 1e-300)

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(8, 6))

        # Scatter plot colored by label
        if label_col in plot_data.columns:
            plot_data[label_col] = plot_data[label_col].astype(str)
            sns.scatterplot(
                data=plot_data, x=log2fc_col, y="neg_log10_pval",
                hue=label_col, alpha=0.7, ax=ax,
            )
            ax.legend(title="Label", bbox_to_anchor=(1.02, 1), loc="upper left")
        else:
            ax.scatter(plot_data[log2fc_col], plot_data["neg_log10_pval"], alpha=0.7)

        # Reference lines
        ax.axhline(-np.log10(pval_threshold), color="red", linestyle="--", alpha=0.5)
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)

        # Label top N guides
        top_n = plot_data.nsmallest(n_label, pvalue_col)
        for _, row in top_n.iterrows():
            ax.text(
                row[log2fc_col], row["neg_log10_pval"],
                row["guide_id"], fontsize=8, alpha=0.8,
            )

        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10(p-value)")
        ax.set_title(f"Intended Target Knockdown (N={len(plot_data)} guides)")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_volcano.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_log2fc_distribution(
    intended_results: pd.DataFrame,
    outdir: str,
    prefix: str = "intended_target",
    log2fc_col: str = "log2_fc",
    label_col: str = "label",
) -> None:
    """
    Plot distribution of log2 fold changes for intended targets.
    """
    plot_data = intended_results.dropna(subset=[log2fc_col]).copy()
    if len(plot_data) == 0:
        logger.warning("No valid data for log2FC distribution plot")
        return

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(10, 5))

        if label_col in plot_data.columns:
            plot_data[label_col] = plot_data[label_col].astype(str)
            sns.histplot(
                data=plot_data, x=log2fc_col, hue=label_col,
                bins=50, ax=ax, element="step", stat="count",
            )
        else:
            sns.histplot(plot_data[log2fc_col], bins=50, ax=ax)

        # Reference lines
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)
        ax.axvline(np.log2(0.4), color="red", linestyle="--", alpha=0.5, label="60% KD")

        median_fc = plot_data[log2fc_col].median()
        ax.axvline(median_fc, color="blue", linestyle="--", alpha=0.5)
        ax.text(0.95, 0.95, f"Median: {median_fc:.2f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=12)

        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of Intended Target Log2 Fold Changes")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_log2fc_distribution.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def run_intended_target_qc(
    input_path: str,
    outdir: str,
    guide_mod_key: str = "guide",
    results_key: str = "trans_per_guide_results",
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
    fc_threshold: float = 0.4,
    pval_threshold: float = 0.05,
    non_targeting_label: str = "non_targeting",
    prefix: str = "intended_target",
) -> None:
    """
    Run intended target inference QC.

    Parameters
    ----------
    input_path : str
        Path to MuData (.h5mu) file.
    outdir : str
        Output directory.
    guide_mod_key : str
        Modality key for guide data.
    results_key : str
        Key in mdata.uns for trans results (default: 'trans_per_guide_results').
        We use trans results because they include all guide-gene pairs, including
        non-targeting guides which are needed for AUROC/AUPRC evaluation.
    log2fc_col : str
        Column name for log2 fold change in results.
    pvalue_col : str
        Column name for p-value in results.
    fc_threshold : float
        Fold change threshold for strong knockdown (default 0.4 = 60% knockdown).
    pval_threshold : float
        P-value threshold for significance.
    non_targeting_label : str
        Value in guide.var["label"] identifying non-targeting guides.
    prefix : str
        Prefix for output filenames.
    """
    os.makedirs(outdir, exist_ok=True)

    # Load data
    logger.info(f"Loading MuData from {input_path}")
    mdata = mudata.read_h5mu(input_path)

    if guide_mod_key not in mdata.mod:
        raise ValueError(f"Modality '{guide_mod_key}' not found in MuData")

    guide = mdata.mod[guide_mod_key]
    logger.info(f"Guide modality: {guide.n_obs} cells, {guide.n_vars} guides")

    # Resolve results key and columns
    resolved_key = resolve_results_key(mdata, results_key)

    # Load inference results (includes all guide-gene pairs)
    trans_results = load_inference_results(mdata, results_key=resolved_key)
    logger.info(f"Loaded {len(trans_results)} test results from '{resolved_key}'")

    # Resolve column names if needed
    log2fc_col, pvalue_col = resolve_metric_columns(
        trans_results, log2fc_col=log2fc_col, pvalue_col=pvalue_col
    )

    # Filter to intended targets (targeting guides only)
    intended = filter_to_intended_targets(
        trans_results, guide.var,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col,
    )
    logger.info(f"Filtered to {len(intended)} intended target tests")

    # Compute knockdown metrics
    metrics = compute_knockdown_metrics(
        intended,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col,
        fc_threshold=fc_threshold, pval_threshold=pval_threshold,
    )

    # -------------------------------------------------------------------------
    # Compute AUROC/AUPRC: targeting vs non-targeting discrimination
    # -------------------------------------------------------------------------
    eval_table = build_evaluation_table(
        trans_results, guide.var,
        pvalue_col=pvalue_col,
        non_targeting_label=non_targeting_label,
    )

    if len(eval_table) > 0:
        auroc, auprc, curves = compute_auroc_auprc(eval_table)
        metrics["auroc"] = auroc
        metrics["auprc"] = auprc
        metrics["n_eval_positives"] = int((eval_table["direct_target"] == 1).sum())
        metrics["n_eval_negatives"] = int((eval_table["direct_target"] == 0).sum())
        logger.info(f"AUROC: {auroc:.3f}, AUPRC: {auprc:.3f}")

        # Plot ROC/PR curves
        plot_roc_pr_curves(curves, auroc, auprc, outdir, prefix=prefix)
    else:
        metrics["auroc"] = np.nan
        metrics["auprc"] = np.nan
        metrics["n_eval_positives"] = 0
        metrics["n_eval_negatives"] = 0
        logger.warning("Could not compute AUROC/AUPRC: insufficient data")

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_path = os.path.join(outdir, f"{prefix}_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info(f"Saved metrics to {metrics_path}")

    # Save full intended target results
    results_path = os.path.join(outdir, f"{prefix}_results.tsv")
    output_cols = ["guide_id", "gene_name", "label", log2fc_col, pvalue_col]
    output_cols = [c for c in output_cols if c in intended.columns]
    intended[output_cols].to_csv(results_path, sep="\t", index=False)
    logger.info(f"Saved results to {results_path}")

    # Generate plots
    plot_volcano(
        intended, outdir, prefix=prefix,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col,
        pval_threshold=pval_threshold,
    )
    plot_log2fc_distribution(
        intended, outdir, prefix=prefix,
        log2fc_col=log2fc_col,
    )

    # Log summary
    logger.info("=" * 50)
    logger.info("Intended Target QC Summary:")
    logger.info(f"  Guides tested: {metrics['n_guides_tested']}")
    logger.info(f"  Strong knockdowns (>={int((1-fc_threshold)*100)}%): "
                f"{metrics['n_strong_knockdowns']} ({metrics['frac_strong_knockdowns']*100:.1f}%)")
    logger.info(f"  Significant (p<{pval_threshold}): "
                f"{metrics['n_significant']} ({metrics['frac_significant']*100:.1f}%)")
    logger.info(f"  Median log2FC: {metrics['median_log2fc']:.3f}")
    if not np.isnan(metrics["auroc"]):
        logger.info(f"  AUROC: {metrics['auroc']:.3f}")
        logger.info(f"  AUPRC: {metrics['auprc']:.3f}")
    logger.info("=" * 50)

    logger.info("Intended target QC complete.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Intended target inference QC metrics and visualizations."
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to MuData (.h5mu) file."
    )
    parser.add_argument(
        "--outdir", "-o", required=True,
        help="Output directory."
    )
    parser.add_argument(
        "--guide-mod-key", default="guide",
        help="Modality key for guide data (default: 'guide')."
    )
    parser.add_argument(
        "--results-key", default="trans_per_guide_results",
        help="Key in mdata.uns for results (default: 'trans_per_guide_results'). "
             "Use 'auto' to pick the best available key."
    )
    parser.add_argument(
        "--log2fc-col", default="log2_fc",
        help="Column name for log2 fold change (default: 'log2_fc'). Use 'auto' to detect."
    )
    parser.add_argument(
        "--pvalue-col", default="p_value",
        help="Column name for p-value (default: 'p_value'). Use 'auto' to detect."
    )
    parser.add_argument(
        "--fc-threshold", type=float, default=0.4,
        help="Fold change threshold for strong knockdown (default: 0.4 = 60%% KD)."
    )
    parser.add_argument(
        "--pval-threshold", type=float, default=0.05,
        help="P-value threshold for significance (default: 0.05)."
    )
    parser.add_argument(
        "--non-targeting-label", default="non_targeting",
        help="Value in guide.var['label'] for non-targeting guides (default: 'non_targeting')."
    )
    parser.add_argument(
        "--prefix", default="intended_target",
        help="Prefix for output filenames (default: 'intended_target')."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_intended_target_qc(
        input_path=args.input,
        outdir=args.outdir,
        guide_mod_key=args.guide_mod_key,
        results_key=args.results_key,
        log2fc_col=args.log2fc_col,
        pvalue_col=args.pvalue_col,
        fc_threshold=args.fc_threshold,
        pval_threshold=args.pval_threshold,
        non_targeting_label=args.non_targeting_label,
        prefix=args.prefix,
    )
