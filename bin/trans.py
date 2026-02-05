#!/usr/bin/env python3
"""
Trans-regulatory inference QC metrics and visualizations.

This script evaluates trans-regulatory effects detected across the genome,
measuring how many genes each perturbation significantly affects beyond its
direct target. Optionally evaluates recovery of validated trans regulatory
relationships using AUROC/AUPRC.

Data structures used:
  - mdata.uns["trans_per_guide_results"]: DataFrame with trans test results
      Columns: guide_id, gene_id, log2_fc, p_value
  - guide.var: Per-guide metadata
      Key columns: guide_id, intended_target_name, gene_name, label

Processing steps:
  1. Load trans_per_guide_results from mdata.uns
  2. Compute per-guide trans metrics (n_significant, median_log2fc, etc.)
  3. Compute overall summary metrics by guide type
  4. If validated_links provided:
     a. Build evaluation table (positives = validated links, negatives = NT guides)
     b. Compute AUROC/AUPRC
     c. Plot ROC/PR curves
  5. Generate volcano plot with subsampled background and highlighted validated links
  6. Generate per-guide distribution plots

Metrics computed:
  Per-guide:
    - n_tests: Total genes tested
    - n_significant_trans: Genes with p_value < threshold (excluding intended target)
    - n_upregulated/n_downregulated: Direction of significant effects
    - median_log2fc, mean_log2fc, max_abs_log2fc
    - frac_significant

  Overall:
    - n_guides_tested, n_targeting_guides, n_non_targeting_guides
    - median_significant_per_guide, total_significant_tests
    - auroc, auprc (if validated links provided)

Outputs:
  - Summary metrics TSV (includes AUROC/AUPRC if validated links)
  - Per-guide summary TSV
  - Volcano plot (subsampled background + highlighted validated links)
  - ROC and Precision-Recall curves (if validated links)
  - Per-guide distribution plots
"""

import argparse
import logging
import os
import sys
from typing import Optional, Tuple, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import mudata
import numpy as np
import pandas as pd
from mudata import MuData
from sklearn.metrics import precision_recall_curve, roc_curve, auc
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

sns.set_style("white")


# ---------------------------------------------------------------------
# Utility Functions
# ---------------------------------------------------------------------
def build_symbol_to_id_map(guide_var: pd.DataFrame) -> pd.Series:
    """
    Build a mapping from gene_name (symbol) to intended_target_name (gene_id).

    Parameters
    ----------
    guide_var : pd.DataFrame
        guide.var with columns "intended_target_name" and "gene_name".

    Returns
    -------
    pd.Series
        Index = gene_name (symbol), values = intended_target_name (gene_id).
    """
    if "intended_target_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'intended_target_name' column")
    if "gene_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'gene_name' column")

    id_map = guide_var.set_index("gene_name")["intended_target_name"]
    id_map = id_map[~id_map.index.duplicated(keep="first")]
    return id_map


def build_id_to_symbol_map(guide_var: pd.DataFrame) -> pd.Series:
    """
    Build a mapping from intended_target_name (gene_id) to gene_name (symbol).

    Parameters
    ----------
    guide_var : pd.DataFrame
        guide.var with columns "intended_target_name" and "gene_name".

    Returns
    -------
    pd.Series
        Index = intended_target_name (gene_id), values = gene_name (symbol).
    """
    if "intended_target_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'intended_target_name' column")
    if "gene_name" not in guide_var.columns:
        raise ValueError("guide.var must have 'gene_name' column")

    id_map = guide_var.set_index("intended_target_name")["gene_name"]
    id_map = id_map[~id_map.index.duplicated(keep="first")]
    return id_map


def load_validated_links(
    validated_links_path: str,
    guide_var: pd.DataFrame,
) -> pd.DataFrame:
    """
    Load and validate trans regulatory links from a TSV file.

    The input TSV has two columns:
      - guide_target: Gene symbol of the perturbed TF (e.g., "TP53")
      - gene: Gene symbol of the expected trans target (e.g., "BEX3")

    Processing steps:
      1. Load TSV with columns [guide_target, gene]
      2. Validate that guide_target exists in guide.var["gene_name"]
      3. Map gene symbols to gene_id (if possible via guide.var mapping)

    Parameters
    ----------
    validated_links_path : str
        Path to TSV file with validated trans links.
    guide_var : pd.DataFrame
        guide.var DataFrame for symbol-to-ID mapping.

    Returns
    -------
    pd.DataFrame
        Validated links with columns: guide_target, gene_symbol.
    """
    links = pd.read_csv(validated_links_path, sep="\t")

    # Validate expected columns
    if "guide_target" not in links.columns or "gene" not in links.columns:
        raise ValueError(
            f"Validated links TSV must have 'guide_target' and 'gene' columns. "
            f"Found: {links.columns.tolist()}"
        )

    # Rename for clarity
    links = links.rename(columns={"gene": "gene_symbol"})

    # Get valid guide targets (TFs that are actually targeted)
    valid_targets = set(guide_var["gene_name"].dropna().unique())

    # Filter to links where guide_target exists in our data
    valid_links = links[links["guide_target"].isin(valid_targets)].copy()

    n_dropped = len(links) - len(valid_links)
    if n_dropped > 0:
        logger.warning(
            f"Dropped {n_dropped} validated links: guide_target not in guide.var['gene_name']"
        )

    logger.info(f"Loaded {len(valid_links)} validated trans links")

    return valid_links


# ---------------------------------------------------------------------
# Multiple Testing Correction
# ---------------------------------------------------------------------
def apply_fdr_correction(
    trans_results: pd.DataFrame,
    pvalue_col: str = "p_value",
    method: str = "fdr_bh",
    alpha: float = 0.05,
) -> pd.DataFrame:
    """
    Apply multiple testing correction to trans results.

    With millions of trans tests, raw p-values will produce many false positives.
    This function applies FDR correction (Benjamini-Hochberg by default) to
    compute adjusted p-values (q-values).

    Parameters
    ----------
    trans_results : pd.DataFrame
        Trans results with p_value column.
    pvalue_col : str
        Column containing raw p-values.
    method : str
        Correction method. Options:
        - "fdr_bh": Benjamini-Hochberg FDR (default, recommended)
        - "bonferroni": Bonferroni correction (very conservative)
        - "none": No correction (use raw p-values)
    alpha : float
        Family-wise error rate for correction.

    Returns
    -------
    pd.DataFrame
        Trans results with added "p_value_adj" column.
    """
    results = trans_results.copy()

    if method == "none":
        results["p_value_adj"] = results[pvalue_col]
        logger.info("No multiple testing correction applied")
        return results

    # Get p-values, handling NaN
    pvals = results[pvalue_col].values
    valid_mask = ~np.isnan(pvals)

    if valid_mask.sum() == 0:
        results["p_value_adj"] = np.nan
        return results

    # Apply correction only to valid p-values
    valid_pvals = pvals[valid_mask]

    logger.info(f"Applying {method} correction to {len(valid_pvals):,} p-values...")

    # multipletests returns: reject, pvals_corrected, alphacSidak, alphacBonf
    _, pvals_adj, _, _ = multipletests(valid_pvals, alpha=alpha, method=method)

    # Put adjusted p-values back
    p_adj = np.full(len(pvals), np.nan)
    p_adj[valid_mask] = pvals_adj
    results["p_value_adj"] = p_adj

    n_sig_raw = (pvals[valid_mask] < alpha).sum()
    n_sig_adj = (pvals_adj < alpha).sum()
    logger.info(
        f"Significant tests: {n_sig_raw:,} (raw p<{alpha}) -> {n_sig_adj:,} (adjusted p<{alpha})"
    )

    return results


# ---------------------------------------------------------------------
# Core Metrics Functions
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

    candidates = [
        "trans_per_guide_results",
        "per_guide_results",
        "cis_per_guide_results",
        "trans_test_results",
        "test_results",
    ]
    for key in candidates:
        if key in available:
            logger.info(f"Using results key '{key}' for trans QC")
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


def compute_per_guide_trans_metrics(
    trans_results: pd.DataFrame,
    guide_var: pd.DataFrame,
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
    pval_threshold: float = 0.05,
    exclude_intended_target: bool = True,
) -> pd.DataFrame:
    """
    Compute per-guide summary statistics for trans effects across the genome.

    For each guide, compute:
      - n_tests: Total number of genes tested
      - n_significant_trans: Number of genes with p_value < threshold
      - n_upregulated: Significant genes with log2fc > 0
      - n_downregulated: Significant genes with log2fc < 0
      - median_log2fc: Median log2FC across all tested genes
      - mean_log2fc: Mean log2FC across all tested genes
      - max_abs_log2fc: Maximum absolute log2FC (strongest effect)
      - frac_significant: Fraction of tests that are significant

    Parameters
    ----------
    trans_results : pd.DataFrame
        Trans results with guide_id, gene_id, log2_fc, p_value.
    guide_var : pd.DataFrame
        guide.var for intended target exclusion and metadata.
    log2fc_col : str
        Column for log2 fold change.
    pvalue_col : str
        Column for p-value.
    pval_threshold : float
        Threshold for significance.
    exclude_intended_target : bool
        If True, exclude intended target from trans metrics (it's a cis effect).

    Returns
    -------
    pd.DataFrame
        Per-guide summary with columns: guide_id, n_tests, n_significant_trans, ...
    """
    results = trans_results.copy()

    # Add gene_name by mapping gene_id
    id_to_symbol = build_id_to_symbol_map(guide_var)
    results["gene_name"] = results["gene_id"].map(id_to_symbol)

    # Extract target_name from guide_id (e.g., "CD81#strong" -> "CD81")
    results["target_name"] = results["guide_id"].str.split("#").str[0]

    if exclude_intended_target:
        # Exclude intended target (cis effect) from trans metrics
        results = results[results["gene_name"] != results["target_name"]].copy()

    # Compute per-guide metrics
    def guide_metrics(df):
        valid = df.dropna(subset=[log2fc_col, pvalue_col])
        n_tests = len(valid)

        if n_tests == 0:
            return pd.Series({
                "n_tests": 0,
                "n_significant_trans": 0,
                "n_upregulated": 0,
                "n_downregulated": 0,
                "median_log2fc": np.nan,
                "mean_log2fc": np.nan,
                "max_abs_log2fc": np.nan,
                "frac_significant": 0.0,
            })

        significant = valid[valid[pvalue_col] < pval_threshold]
        n_sig = len(significant)
        n_up = len(significant[significant[log2fc_col] > 0])
        n_down = len(significant[significant[log2fc_col] < 0])

        return pd.Series({
            "n_tests": n_tests,
            "n_significant_trans": n_sig,
            "n_upregulated": n_up,
            "n_downregulated": n_down,
            "median_log2fc": valid[log2fc_col].median(),
            "mean_log2fc": valid[log2fc_col].mean(),
            "max_abs_log2fc": valid[log2fc_col].abs().max(),
            "frac_significant": n_sig / n_tests if n_tests > 0 else 0.0,
        })

    per_guide = results.groupby("guide_id").apply(guide_metrics).reset_index()

    # Merge with guide metadata
    guide_meta = guide_var.reset_index(drop=True)[["guide_id", "gene_name", "label"]].drop_duplicates()
    per_guide = per_guide.merge(guide_meta, on="guide_id", how="left")

    return per_guide


def compute_overall_trans_metrics(
    per_guide_summary: pd.DataFrame,
    non_targeting_label: str = "non_targeting",
) -> dict:
    """
    Compute overall summary metrics across all guides.

    Parameters
    ----------
    per_guide_summary : pd.DataFrame
        Output from compute_per_guide_trans_metrics.
    non_targeting_label : str
        Label for non-targeting guides.

    Returns
    -------
    dict
        Dictionary of overall metrics.
    """
    targeting = per_guide_summary[per_guide_summary["label"] != non_targeting_label]
    non_targeting = per_guide_summary[per_guide_summary["label"] == non_targeting_label]

    metrics = {
        "n_guides_tested": len(per_guide_summary),
        "n_targeting_guides": len(targeting),
        "n_non_targeting_guides": len(non_targeting),
        "median_significant_per_guide_targeting": targeting["n_significant_trans"].median() if len(targeting) > 0 else np.nan,
        "mean_significant_per_guide_targeting": targeting["n_significant_trans"].mean() if len(targeting) > 0 else np.nan,
        "median_significant_per_guide_nt": non_targeting["n_significant_trans"].median() if len(non_targeting) > 0 else np.nan,
        "mean_significant_per_guide_nt": non_targeting["n_significant_trans"].mean() if len(non_targeting) > 0 else np.nan,
        "total_significant_tests": int(per_guide_summary["n_significant_trans"].sum()),
        "median_genome_log2fc_targeting": targeting["median_log2fc"].median() if len(targeting) > 0 else np.nan,
        "median_genome_log2fc_nt": non_targeting["median_log2fc"].median() if len(non_targeting) > 0 else np.nan,
    }

    return metrics


# ---------------------------------------------------------------------
# AUROC/AUPRC Evaluation
# ---------------------------------------------------------------------
def build_validated_links_evaluation_table(
    trans_results: pd.DataFrame,
    guide_var: pd.DataFrame,
    validated_links: pd.DataFrame,
    pvalue_col: str = "p_value",
    non_targeting_label: str = "non_targeting",
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Build evaluation table for AUROC/AUPRC calculation using validated trans links.

    Positive controls (validated_link=1):
      - Tests where a TF-targeting guide affects a validated trans target gene
      - Match: guide's target TF == validated_links["guide_target"] AND
               gene matches validated_links["gene_symbol"]

    Negative controls (validated_link=0):
      - Non-targeting guides tested against the same validated target genes
      - Provides matched controls: same genes, but no expected effect

    The dataset is balanced by downsampling negative controls.

    Parameters
    ----------
    trans_results : pd.DataFrame
        Trans test results.
    guide_var : pd.DataFrame
        guide.var with label, gene_name columns.
    validated_links : pd.DataFrame
        Output from load_validated_links with columns: guide_target, gene_symbol.
    pvalue_col : str
        Column containing p-values.
    non_targeting_label : str
        Label for non-targeting guides.
    random_state : int
        Random seed for balanced sampling.

    Returns
    -------
    pd.DataFrame
        Evaluation table with columns: guide_id, gene_id, p_value, validated_link.
    """
    # Build ID mappings
    id_to_symbol = build_id_to_symbol_map(guide_var)

    # Prepare results with gene_name and target_name
    results = trans_results.copy()
    results["gene_name"] = results["gene_id"].map(id_to_symbol)
    results["target_name"] = results["guide_id"].str.split("#").str[0]

    # Merge guide metadata to get label
    guide_meta = guide_var.reset_index(drop=True)[["guide_id", "label"]].drop_duplicates()
    results = results.merge(guide_meta, on="guide_id", how="left")

    # Get set of validated target genes
    validated_target_genes = set(validated_links["gene_symbol"].unique())

    # -------------------------------------------------------------------------
    # Positive controls: targeting guide tests validated trans targets
    # -------------------------------------------------------------------------
    # Create a set of (guide_target, gene_symbol) pairs for fast lookup
    validated_pairs = set(zip(validated_links["guide_target"], validated_links["gene_symbol"]))

    # For each row, check if (target_name, gene_name) is in validated_pairs
    results["is_validated_pair"] = results.apply(
        lambda r: (r["target_name"], r["gene_name"]) in validated_pairs
        if pd.notna(r["target_name"]) and pd.notna(r["gene_name"])
        else False,
        axis=1,
    )

    positives = results[
        (results["is_validated_pair"]) &
        (results["label"] != non_targeting_label) &
        (results[pvalue_col].notna())
    ].copy()
    positives["validated_link"] = 1

    # -------------------------------------------------------------------------
    # Negative controls: non-targeting guides tested against validated target genes
    # -------------------------------------------------------------------------
    negatives = results[
        (results["label"] == non_targeting_label) &
        (results["gene_name"].isin(validated_target_genes)) &
        (results[pvalue_col].notna())
    ].copy()
    negatives["validated_link"] = 0

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
    eval_table = eval_table[["guide_id", "gene_id", pvalue_col, "validated_link"]].copy()
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
        Must have columns: p_value, validated_link.

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

    true_labels = eval_table["validated_link"].values
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
# Plotting Functions
# ---------------------------------------------------------------------
def plot_roc_pr_curves(
    curves: dict,
    auroc: float,
    auprc: float,
    outdir: str,
    prefix: str = "trans",
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


def plot_trans_volcano(
    trans_results: pd.DataFrame,
    outdir: str,
    validated_links: Optional[pd.DataFrame] = None,
    guide_var: Optional[pd.DataFrame] = None,
    prefix: str = "trans",
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
    pval_threshold: float = 0.05,
    n_background: int = 10000,
    random_state: int = 42,
) -> None:
    """
    Plot volcano plot of trans effects with subsampled background.

    Plotting approach:
      1. Subsample n_background random trans pairs as grey background points
      2. If validated_links provided, overlay validated link tests in red

    Parameters
    ----------
    trans_results : pd.DataFrame
        All trans test results.
    outdir : str
        Output directory.
    validated_links : pd.DataFrame, optional
        Validated trans links (from load_validated_links).
    guide_var : pd.DataFrame, optional
        Required if validated_links is provided, for target mapping.
    prefix : str
        Filename prefix.
    log2fc_col, pvalue_col : str
        Column names.
    pval_threshold : float
        Significance threshold for reference line.
    n_background : int
        Number of background points to sample.
    random_state : int
        Random seed for reproducibility.
    """
    plot_data = trans_results.dropna(subset=[log2fc_col, pvalue_col]).copy()
    if len(plot_data) == 0:
        logger.warning("No valid data for volcano plot")
        return

    # Compute -log10(p-value)
    plot_data["neg_log10_pval"] = -np.log10(plot_data[pvalue_col] + 1e-300)

    # Subsample background
    if len(plot_data) > n_background:
        background = plot_data.sample(n=n_background, random_state=random_state)
    else:
        background = plot_data

    with sns.plotting_context("talk"):
        fig, ax = plt.subplots(figsize=(10, 7))

        # Plot background points in grey
        ax.scatter(
            background[log2fc_col],
            background["neg_log10_pval"],
            c="lightgrey",
            alpha=0.5,
            s=10,
            label=f"Background (n={len(background)})",
        )

        # Overlay validated links if provided
        if validated_links is not None and guide_var is not None:
            id_to_symbol = build_id_to_symbol_map(guide_var)
            plot_data["gene_name"] = plot_data["gene_id"].map(id_to_symbol)
            plot_data["target_name"] = plot_data["guide_id"].str.split("#").str[0]

            # Find validated pairs
            validated_pairs = set(zip(validated_links["guide_target"], validated_links["gene_symbol"]))

            validated_mask = plot_data.apply(
                lambda r: (r["target_name"], r["gene_name"]) in validated_pairs
                if pd.notna(r["target_name"]) and pd.notna(r["gene_name"])
                else False,
                axis=1,
            )

            validated_data = plot_data[validated_mask]

            if len(validated_data) > 0:
                ax.scatter(
                    validated_data[log2fc_col],
                    validated_data["neg_log10_pval"],
                    c="#e74c3c",
                    alpha=0.8,
                    s=30,
                    label=f"Validated links (n={len(validated_data)})",
                    zorder=10,
                )

        # Reference lines
        ax.axhline(-np.log10(pval_threshold), color="blue", linestyle="--", alpha=0.5, label=f"p={pval_threshold}")
        ax.axvline(0, color="grey", linestyle="--", alpha=0.5)

        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10(p-value)")
        ax.set_title(f"Trans Effects Volcano (Total: {len(trans_results):,} tests)")
        ax.legend(loc="upper right")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_volcano.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


def plot_per_guide_trans_distribution(
    per_guide_summary: pd.DataFrame,
    outdir: str,
    prefix: str = "trans",
    label_col: str = "label",
) -> None:
    """
    Plot distribution of n_significant_trans per guide, colored by guide label.

    Two panels:
      1. Histogram of n_significant_trans, colored by label
      2. Boxplot of n_significant_trans by label category

    Parameters
    ----------
    per_guide_summary : pd.DataFrame
        Per-guide summary from compute_per_guide_trans_metrics.
    outdir : str
        Output directory.
    prefix : str
        Filename prefix.
    label_col : str
        Column for coloring.
    """
    plot_data = per_guide_summary.dropna(subset=["n_significant_trans"]).copy()
    if len(plot_data) == 0:
        logger.warning("No valid data for per-guide distribution plot")
        return

    with sns.plotting_context("talk"):
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # Histogram
        ax = axes[0]
        if label_col in plot_data.columns:
            plot_data[label_col] = plot_data[label_col].astype(str)
            sns.histplot(
                data=plot_data, x="n_significant_trans", hue=label_col,
                bins=50, ax=ax, element="step", stat="count",
            )
        else:
            sns.histplot(plot_data["n_significant_trans"], bins=50, ax=ax)

        ax.set_xlabel("Number of Significant Trans Effects")
        ax.set_ylabel("Count")
        ax.set_title("Distribution of Trans Effects per Guide")

        # Boxplot
        ax = axes[1]
        if label_col in plot_data.columns:
            order = sorted(plot_data[label_col].unique())
            sns.boxplot(
                data=plot_data, x=label_col, y="n_significant_trans",
                order=order, ax=ax,
            )
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        else:
            sns.boxplot(y=plot_data["n_significant_trans"], ax=ax)

        ax.set_xlabel("Guide Label")
        ax.set_ylabel("Number of Significant Trans Effects")
        ax.set_title("Trans Effects by Guide Type")

        plt.tight_layout()
        fname = os.path.join(outdir, f"{prefix}_per_guide_distribution.png")
        plt.savefig(fname, dpi=200)
        plt.close()
        logger.info(f"Saved {fname}")


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------
def run_trans_qc(
    input_path: str,
    outdir: str,
    validated_links_path: Optional[str] = None,
    guide_mod_key: str = "guide",
    results_key: str = "trans_per_guide_results",
    log2fc_col: str = "log2_fc",
    pvalue_col: str = "p_value",
    pval_threshold: float = 0.05,
    fdr_method: str = "fdr_bh",
    non_targeting_label: str = "non_targeting",
    n_background: int = 10000,
    prefix: str = "trans",
) -> None:
    """
    Run trans-regulatory inference QC.

    Parameters
    ----------
    input_path : str
        Path to MuData (.h5mu) file.
    outdir : str
        Output directory.
    validated_links_path : str, optional
        Path to TSV with validated trans links. If None, skip AUROC/AUPRC.
    guide_mod_key : str
        Modality key for guide data.
    results_key : str
        Key in mdata.uns for trans results.
    log2fc_col : str
        Column name for log2 fold change.
    pvalue_col : str
        Column name for p-value (raw).
    pval_threshold : float
        P-value threshold for significance (applied to adjusted p-values).
    fdr_method : str
        Multiple testing correction method: "fdr_bh" (Benjamini-Hochberg),
        "bonferroni", or "none" (no correction).
    non_targeting_label : str
        Value in guide.var["label"] for non-targeting guides.
    n_background : int
        Number of background points for volcano plot.
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

    # Load inference results
    trans_results = load_inference_results(mdata, results_key=resolved_key)
    logger.info(f"Loaded {len(trans_results):,} test results from '{resolved_key}'")

    # Resolve column names if needed
    log2fc_col, pvalue_col = resolve_metric_columns(
        trans_results, log2fc_col=log2fc_col, pvalue_col=pvalue_col
    )

    # Apply multiple testing correction
    trans_results = apply_fdr_correction(
        trans_results, pvalue_col=pvalue_col, method=fdr_method, alpha=pval_threshold
    )
    pvalue_col_adj = "p_value_adj"  # Use adjusted p-values for significance

    # Compute per-guide trans metrics (using adjusted p-values for significance)
    per_guide_summary = compute_per_guide_trans_metrics(
        trans_results, guide.var,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col_adj,
        pval_threshold=pval_threshold,
    )
    logger.info(f"Computed per-guide metrics for {len(per_guide_summary)} guides")

    # Compute overall metrics
    metrics = compute_overall_trans_metrics(per_guide_summary, non_targeting_label=non_targeting_label)
    metrics["fdr_method"] = fdr_method

    # -------------------------------------------------------------------------
    # Optional: AUROC/AUPRC evaluation with validated links
    # -------------------------------------------------------------------------
    validated_links = None
    if validated_links_path:
        logger.info(f"Loading validated trans links from {validated_links_path}")
        validated_links = load_validated_links(validated_links_path, guide.var)

        if len(validated_links) > 0:
            eval_table = build_validated_links_evaluation_table(
                trans_results, guide.var, validated_links,
                pvalue_col=pvalue_col,
                non_targeting_label=non_targeting_label,
            )

            if len(eval_table) > 0:
                auroc, auprc, curves = compute_auroc_auprc(eval_table)
                metrics["auroc"] = auroc
                metrics["auprc"] = auprc
                metrics["n_validated_links"] = len(validated_links)
                metrics["n_eval_positives"] = int((eval_table["validated_link"] == 1).sum())
                metrics["n_eval_negatives"] = int((eval_table["validated_link"] == 0).sum())
                logger.info(f"AUROC: {auroc:.3f}, AUPRC: {auprc:.3f}")

                # Plot ROC/PR curves
                plot_roc_pr_curves(curves, auroc, auprc, outdir, prefix=prefix)

                # Save validated links results
                validated_results = trans_results.copy()
                id_to_symbol = build_id_to_symbol_map(guide.var)
                validated_results["gene_name"] = validated_results["gene_id"].map(id_to_symbol)
                validated_results["target_name"] = validated_results["guide_id"].str.split("#").str[0]

                validated_pairs = set(zip(validated_links["guide_target"], validated_links["gene_symbol"]))
                validated_mask = validated_results.apply(
                    lambda r: (r["target_name"], r["gene_name"]) in validated_pairs
                    if pd.notna(r["target_name"]) and pd.notna(r["gene_name"])
                    else False,
                    axis=1,
                )
                validated_results = validated_results[validated_mask].copy()
                validated_results["significant"] = validated_results[pvalue_col_adj] < pval_threshold

                validated_path = os.path.join(outdir, f"{prefix}_validated_links_results.tsv")
                validated_results.to_csv(validated_path, sep="\t", index=False)
                logger.info(f"Saved validated links results to {validated_path}")
            else:
                logger.warning("Could not compute AUROC/AUPRC: insufficient data")
                metrics["auroc"] = np.nan
                metrics["auprc"] = np.nan
                metrics["n_validated_links"] = len(validated_links)
                metrics["n_eval_positives"] = 0
                metrics["n_eval_negatives"] = 0
    else:
        metrics["auroc"] = np.nan
        metrics["auprc"] = np.nan
        metrics["n_validated_links"] = 0
        metrics["n_eval_positives"] = 0
        metrics["n_eval_negatives"] = 0

    # Save metrics
    metrics_df = pd.DataFrame([metrics])
    metrics_path = os.path.join(outdir, f"{prefix}_metrics.tsv")
    metrics_df.to_csv(metrics_path, sep="\t", index=False)
    logger.info(f"Saved metrics to {metrics_path}")

    # Save per-guide summary
    per_guide_path = os.path.join(outdir, f"{prefix}_per_guide_summary.tsv")
    per_guide_summary.to_csv(per_guide_path, sep="\t", index=False)
    logger.info(f"Saved per-guide summary to {per_guide_path}")

    # Save full trans results (with annotations)
    # Add clear columns for guide target and tested gene (both IDs and symbols)
    trans_results_annotated = trans_results.copy()

    # Get gene modality for tested gene symbol lookup
    gene = mdata.mod["gene"]
    gene_id_to_symbol = gene.var["symbol"].to_dict()

    # Guide target: extract symbol from guide_id, get ID from guide.var
    trans_results_annotated["guide_target_symbol"] = trans_results_annotated["guide_id"].str.split("#").str[0]

    # Build symbol -> ID mapping from guide.var for guide targets
    guide_symbol_to_id = guide.var.set_index("gene_name")["intended_target_name"].to_dict()
    trans_results_annotated["guide_target_id"] = trans_results_annotated["guide_target_symbol"].map(guide_symbol_to_id)

    # Tested gene: rename gene_id for clarity, add symbol from gene.var
    trans_results_annotated = trans_results_annotated.rename(columns={"gene_id": "tested_gene_id"})
    trans_results_annotated["tested_gene_symbol"] = trans_results_annotated["tested_gene_id"].map(gene_id_to_symbol)

    # Merge guide label
    guide_meta = guide.var.reset_index(drop=True)[["guide_id", "label"]].drop_duplicates()
    trans_results_annotated = trans_results_annotated.merge(guide_meta, on="guide_id", how="left")

    # Add significant column
    trans_results_annotated["significant"] = trans_results_annotated[pvalue_col_adj] < pval_threshold

    # Define output columns with clear naming
    output_cols = [
        "guide_id",
        "guide_target_id", "guide_target_symbol",
        "tested_gene_id", "tested_gene_symbol",
        "label",
        log2fc_col, pvalue_col, pvalue_col_adj, "significant"
    ]
    output_cols = [c for c in output_cols if c in trans_results_annotated.columns]

    # Save all results
    all_results_path = os.path.join(outdir, f"{prefix}_results.tsv")
    trans_results_annotated[output_cols].to_csv(all_results_path, sep="\t", index=False)
    logger.info(f"Saved {len(trans_results_annotated):,} trans results to {all_results_path}")

    # Save only significant results (smaller file for quick access)
    sig_results = trans_results_annotated[trans_results_annotated["significant"]].copy()
    sig_results_path = os.path.join(outdir, f"{prefix}_significant_results.tsv")
    sig_results[output_cols].to_csv(sig_results_path, sep="\t", index=False)
    logger.info(f"Saved {len(sig_results):,} significant trans results to {sig_results_path}")

    # Generate plots (use adjusted p-values for volcano)
    plot_trans_volcano(
        trans_results, outdir,
        validated_links=validated_links,
        guide_var=guide.var if validated_links is not None else None,
        prefix=prefix,
        log2fc_col=log2fc_col, pvalue_col=pvalue_col_adj,
        pval_threshold=pval_threshold,
        n_background=n_background,
    )

    plot_per_guide_trans_distribution(
        per_guide_summary, outdir, prefix=prefix,
    )

    # Log summary
    logger.info("=" * 50)
    logger.info("Trans-Regulatory QC Summary:")
    logger.info(f"  Guides tested: {metrics['n_guides_tested']}")
    logger.info(f"  Targeting guides: {metrics['n_targeting_guides']}")
    logger.info(f"  Non-targeting guides: {metrics['n_non_targeting_guides']}")
    logger.info(f"  Median sig trans per targeting guide: {metrics['median_significant_per_guide_targeting']:.1f}")
    logger.info(f"  Median sig trans per NT guide: {metrics['median_significant_per_guide_nt']:.1f}")
    logger.info(f"  Total significant trans tests: {metrics['total_significant_tests']:,}")
    if not np.isnan(metrics["auroc"]):
        logger.info(f"  Validated links: {metrics['n_validated_links']}")
        logger.info(f"  AUROC: {metrics['auroc']:.3f}")
        logger.info(f"  AUPRC: {metrics['auprc']:.3f}")
    logger.info("=" * 50)

    logger.info("Trans-regulatory QC complete.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Trans-regulatory inference QC metrics and visualizations."
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
        "--validated-links",
        help="Path to TSV with validated trans links (optional). "
             "Format: two columns 'guide_target' and 'gene' with gene symbols."
    )
    parser.add_argument(
        "--guide-mod-key", default="guide",
        help="Modality key for guide data (default: 'guide')."
    )
    parser.add_argument(
        "--results-key", default="trans_per_guide_results",
        help="Key in mdata.uns for trans results (default: 'trans_per_guide_results'). Use 'auto' to detect."
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
        "--pval-threshold", type=float, default=0.05,
        help="P-value threshold for significance (default: 0.05). Applied to adjusted p-values."
    )
    parser.add_argument(
        "--fdr-method", default="fdr_bh",
        choices=["fdr_bh", "bonferroni", "none"],
        help="Multiple testing correction method (default: 'fdr_bh'). "
             "Options: 'fdr_bh' (Benjamini-Hochberg FDR), 'bonferroni', 'none'."
    )
    parser.add_argument(
        "--non-targeting-label", default="non_targeting",
        help="Value in guide.var['label'] for non-targeting guides (default: 'non_targeting')."
    )
    parser.add_argument(
        "--n-background", type=int, default=10000,
        help="Number of background points for volcano plot (default: 10000)."
    )
    parser.add_argument(
        "--prefix", default="trans",
        help="Prefix for output filenames (default: 'trans')."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_trans_qc(
        input_path=args.input,
        outdir=args.outdir,
        validated_links_path=args.validated_links,
        guide_mod_key=args.guide_mod_key,
        results_key=args.results_key,
        log2fc_col=args.log2fc_col,
        pvalue_col=args.pvalue_col,
        pval_threshold=args.pval_threshold,
        fdr_method=args.fdr_method,
        non_targeting_label=args.non_targeting_label,
        n_background=args.n_background,
        prefix=args.prefix,
    )
