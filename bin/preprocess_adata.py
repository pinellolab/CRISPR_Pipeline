#!/usr/bin/env python3
"""Apply RNA cell QC to one measurement set before concatenation."""

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata_concat import apply_batch_suffix, get_barcode_key
from count_matrix_utils import normalize_sparse_index_dtypes, to_sparse_counts
from scipy.interpolate import UnivariateSpline


def get_vectors(x, y):
    smooth_spline = UnivariateSpline(x, y, s=len(x))
    second_deriv = smooth_spline.derivative(n=2)(x)
    ten_percent = max(1, round(len(x) * 0.1))
    mid = second_deriv[ten_percent:-ten_percent] if len(second_deriv) > 2 * ten_percent else second_deriv
    if np.all(mid >= 0) or np.all(mid <= 0):
        return x, y
    minimum = int(np.argmin(second_deriv))
    left = np.where(second_deriv[: minimum + 1] >= 0)[0]
    right = np.where(second_deriv[minimum:] >= 0)[0]
    if not len(left) or not len(right):
        return x, y
    start, end = left[-1], minimum + right[0]
    return (x, y) if start >= end else (x[start : end + 1], y[start : end + 1])


def elbow_knee_finder(x, y, mode="basic"):
    if mode == "advanced":
        if len(np.unique(x)) < 4:
            return None
        x, y = get_vectors(x, y)
    if len(x) == 0 or len(y) == 0 or x[0] == x[-1]:
        return None
    slope = (y[-1] - y[0]) / (x[-1] - x[0])
    intercept = y[0] - slope * x[0]
    distances = np.abs(slope * x - y + intercept) / np.sqrt(slope**2 + 1)
    index = int(np.argmax(distances))
    return np.array([x[index], y[index]])


def get_elbow_knee_points(x, y):
    point_1 = elbow_knee_finder(x, y, mode="basic")
    point_2 = None
    if point_1 is not None:
        end = max(1, min(len(x), int(round(point_1[0]))))
        point_2 = elbow_knee_finder(x[:end], y[:end], mode="advanced")
    return point_1, point_2


def safe_label(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "measurement_set"


def batch_name(mapping_dir):
    match = re.match(r"(.+)_ks_transcripts_out$", Path(mapping_dir).name)
    return match.group(1) if match else Path(mapping_dir).name


def median_mad(values):
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if not finite.size:
        return np.nan, np.nan
    median = float(np.median(finite))
    return median, float(np.median(np.abs(finite - median)))


def mad_limits(values, n_mads, upper_only=False):
    median, mad = median_mad(values)
    if n_mads <= 0 or not np.isfinite(mad) or mad == 0:
        return median, mad, -np.inf, np.inf
    lower = -np.inf if upper_only else median - n_mads * mad
    return median, mad, lower, median + n_mads * mad


def plot_barcode_rank(
    knee_df, point_1, point_2, selected_rank, batch, outpath,
    barcode_filter="knee", selected_threshold=np.nan,
):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(knee_df["rank"], knee_df["sum_log"], linewidth=1.2, color="#2563eb")
    if selected_rank is not None:
        ax.axvline(selected_rank, color="#dc2626", linestyle="--", label="Applied knee")
        ax.axvspan(1, selected_rank, color="#22c55e", alpha=0.10, label="Retained barcodes")
    if barcode_filter == "none":
        note = f"QC_barcode_filter = none\nNo knee filter applied\n{len(knee_df):,} barcodes retained"
    elif selected_rank is None:
        note = f"QC_barcode_filter = {barcode_filter}\nNo valid knee found; filter skipped"
    else:
        note = (
            f"QC_barcode_filter = {barcode_filter}\n"
            f"selected rank = {selected_rank:,}\nRNA UMI threshold ≥ {selected_threshold:,.0f}"
        )
    ax.text(
        0.98, 0.97, note, transform=ax.transAxes, ha="right", va="top", fontsize=9,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.95},
    )
    ax.set(xlabel="Barcode rank", ylabel="Log1p RNA UMI counts", title=f"{batch}: applied barcode-rank knee")
    if selected_rank is not None:
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def plot_filter_steps(snapshots, batch, outpath, knee_context):
    """Plot the applied knee first, then every downstream cell filter."""
    row_count = len(snapshots) + 1
    fig, axes = plt.subplots(row_count, 3, figsize=(17, 3.15 * row_count), squeeze=False)

    knee_df = knee_context["knee_df"]
    selected_rank = knee_context["selected_rank"]
    selected_threshold = knee_context["selected_threshold"]
    cells_before = knee_context["cells_before"]
    cells_after = knee_context["cells_after"]
    applied = selected_rank is not None
    tie_note = "\n(includes UMI-threshold ties)" if applied and cells_after != selected_rank else ""

    full_ax = axes[0, 0]
    full_ax.plot(knee_df["rank"], knee_df["sum_log"], color="#2563eb", linewidth=1.2)
    if applied:
        full_ax.axvline(selected_rank, color="#dc2626", linestyle="--", label="Applied knee")
        full_ax.axvspan(1, cells_after, color="#22c55e", alpha=0.10, label="Retained")
        full_ax.legend(frameon=False, fontsize=8)
    full_ax.set(
        title=f"Before knee: {cells_before:,} barcodes",
        xlabel="Barcode rank",
        ylabel="Log1p RNA UMIs",
    )

    retained_ax = axes[0, 1]
    retained_end = cells_after if applied else len(knee_df)
    retained_curve = knee_df.iloc[:retained_end]
    retained_ax.plot(retained_curve["rank"], retained_curve["sum_log"], color="#16a34a", linewidth=1.2)
    if applied:
        retained_ax.axvline(selected_rank, color="#dc2626", linestyle="--", label="Applied knee")
        retained_ax.legend(frameon=False, fontsize=8)
    retained_ax.set(
        title=f"After knee: {cells_after:,} cells retained{tie_note}",
        xlabel="Retained barcode rank",
        ylabel="Log1p RNA UMIs",
    )

    summary_ax = axes[0, 2]
    summary_ax.axis("off")
    status = "APPLIED" if applied else "DISABLED / SKIPPED"
    threshold = (
        f"rank ≤ {selected_rank:,}\nRNA UMIs ≥ {selected_threshold:,.0f}"
        if applied else "No applied knee"
    )
    summary_ax.text(
        0.5, 0.78, f"{cells_before:,} barcodes before", ha="center", va="center",
        fontsize=12, fontweight="bold",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#dbeafe", "edgecolor": "#2563eb"},
    )
    summary_ax.annotate("", xy=(0.5, 0.56), xytext=(0.5, 0.69), arrowprops={"arrowstyle": "->", "color": "#64748b"})
    summary_ax.text(
        0.5, 0.47, f"QC_barcode_filter = {knee_context['barcode_filter']}\n{threshold}\n{status}",
        ha="center", va="center", fontsize=9.5,
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f8fafc", "edgecolor": "#94a3b8"},
    )
    summary_ax.annotate("", xy=(0.5, 0.24), xytext=(0.5, 0.34), arrowprops={"arrowstyle": "->", "color": "#64748b"})
    summary_ax.text(
        0.5, 0.14,
        f"{cells_after:,} cells after{tie_note}\n{cells_before - cells_after:,} removed",
        ha="center", va="center", fontsize=11, fontweight="bold",
        bbox={"boxstyle": "round,pad=0.45", "facecolor": "#dcfce7", "edgecolor": "#16a34a"},
    )

    for row_index, snapshot in enumerate(snapshots):
        row_index += 1
        before = np.asarray(snapshot["values_before"], dtype=float)
        after = np.asarray(snapshot["values_after"], dtype=float)
        before = before[np.isfinite(before)]
        after = after[np.isfinite(after)]
        finite = before if before.size else np.array([0.0])
        low, high = np.nanpercentile(finite, [0.5, 99.5]) if finite.size > 1 else (finite[0] - 0.5, finite[0] + 0.5)
        if not np.isfinite(low) or not np.isfinite(high) or low == high:
            low, high = float(np.nanmin(finite)) - 0.5, float(np.nanmax(finite)) + 0.5
        span = max(high - low, 1e-6)
        display_xlim = (low - 0.08 * span, high + 0.08 * span)
        bins = np.linspace(low, high, 51)
        for column, (values, color, state) in enumerate(
            ((before, "#60a5fa", "Before"), (after, "#34d399", "After"))
        ):
            ax = axes[row_index, column]
            ax.hist(values, bins=bins, color=color, edgecolor="white")
            # Keep the first panel free for its filter/removal annotation; the
            # identical cutoff legend is shown on the after-filter panel.
            add_filter_bounds(ax, snapshot, display_xlim=display_xlim, show_legend=(column == 1))
            ax.set_xlim(display_xlim)
            count_key = "cells_before" if state == "Before" else "cells_after"
            ax.set_title(f"{state}: {snapshot[count_key]:,} cells")
            ax.set_xlabel(snapshot["metric_label"])
            ax.set_ylabel("Cells")
        box_ax = axes[row_index, 2]
        box_values = [before if before.size else np.array([np.nan]), after if after.size else np.array([np.nan])]
        boxes = box_ax.boxplot(
            box_values, vert=False, tick_labels=["Before", "After"], showfliers=False,
            patch_artist=True, medianprops={"color": "#111827", "linewidth": 1.5},
        )
        for patch, color in zip(boxes["boxes"], ("#93c5fd", "#6ee7b7")):
            patch.set_facecolor(color)
        add_filter_bounds(box_ax, snapshot, display_xlim=display_xlim, show_legend=True)
        box_ax.set_xlim(display_xlim)
        box_ax.set_title("Before/after distribution")
        box_ax.set_xlabel(snapshot["metric_label"])
        status = "applied" if snapshot["applied"] else "disabled/skipped"
        axes[row_index, 0].text(
            0.01, 0.96,
            f"{snapshot['filter_label']} ({status})\nRemoved: {snapshot['removed']:,}",
            transform=axes[row_index, 0].transAxes, ha="left", va="top", fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#cbd5e1", "alpha": 0.92},
        )
    fig.suptitle(f"{batch}: sequential RNA cell-filter impact", fontsize=15, y=0.999)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, facecolor="white")
    plt.close(fig)


def add_filter_bounds(ax, snapshot, show_legend=False, display_xlim=None):
    """Draw filter limits without allowing distant bounds to flatten the data."""
    if not snapshot["applied"]:
        return
    colors = {"MAD": "#dc2626", "Scrublet": "#7c3aed", "fixed": "#059669"}
    color = colors.get(snapshot["bound_kind"], "#059669")
    linestyle = "--" if snapshot["bound_kind"] in {"MAD", "Scrublet"} else ":"
    for bound, label in (
        (snapshot["lower"], snapshot.get("lower_label")),
        (snapshot["upper"], snapshot.get("upper_label")),
    ):
        if np.isfinite(bound):
            resolved_label = f"{label} = {bound:.3g}" if label else f"Threshold = {bound:.3g}"
            outside_view = display_xlim is not None and not (display_xlim[0] <= bound <= display_xlim[1])
            if outside_view:
                resolved_label += " (outside view)"
                # An empty artist preserves the exact cutoff in the legend without
                # stretching the axis until the observed histogram disappears.
                ax.plot([], [], color=color, linestyle=linestyle, linewidth=1.5, label=resolved_label)
            else:
                ax.axvline(bound, color=color, linestyle=linestyle, linewidth=1.5, label=resolved_label)
    if show_legend and (snapshot.get("lower_label") or snapshot.get("upper_label")):
        ax.legend(frameon=False, fontsize=8, loc="best")


def plot_filter_flow(flow, batch, outpath):
    """Render one explicit before → filter → after row per QC step."""
    fig, axes = plt.subplots(len(flow), 1, figsize=(16, max(4, 2.25 * len(flow))), squeeze=False)
    for row_index, (_, step) in enumerate(flow.iterrows()):
        ax = axes[row_index, 0]
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        status = "APPLIED" if bool(step["applied"]) else "DISABLED / SKIPPED"
        before = int(step["cells_before"])
        after = int(step["cells_after"])
        removed = int(step["cells_removed"])
        after_color = "#dcfce7" if removed == 0 else "#fef3c7"
        ax.text(
            0.12, 0.5, f"{before:,}\ncells before", ha="center", va="center",
            fontsize=11, fontweight="bold",
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "#dbeafe", "edgecolor": "#2563eb"},
        )
        ax.annotate("", xy=(0.31, 0.5), xytext=(0.21, 0.5), arrowprops={"arrowstyle": "->", "color": "#64748b"})
        ax.text(
            0.5, 0.5,
            f"Step {int(step['step_order'])}: {step['filter_label']}\n{step['threshold']}\n{status}",
            ha="center", va="center", fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f8fafc", "edgecolor": "#94a3b8"},
        )
        ax.annotate("", xy=(0.79, 0.5), xytext=(0.69, 0.5), arrowprops={"arrowstyle": "->", "color": "#64748b"})
        ax.text(
            0.88, 0.5,
            f"{after:,}\ncells after\n{removed:,} removed ({float(step['removed_percent']):.2f}%)",
            ha="center", va="center", fontsize=10.5, fontweight="bold",
            bbox={"boxstyle": "round,pad=0.5", "facecolor": after_color, "edgecolor": "#16a34a"},
        )
    fig.suptitle(f"{batch}: RNA QC filtering flow (pipeline order)", fontsize=15)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_qc_distributions(obs, limits, fixed_limits, mad_counts, batch, outpath):
    specs = [
        ("log1p_total_counts", "Log1p total RNA UMIs", limits["total_counts"]),
        ("log1p_n_genes_by_counts", "Log1p detected genes", limits["n_genes"]),
        ("pct_counts_mt", "Mitochondrial counts (%)", (np.nan, np.nan, -np.inf, np.inf)),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, (column, label, (_median, _mad, lower, upper)) in zip(axes, specs):
        values = pd.to_numeric(obs[column], errors="coerce").dropna()
        ax.hist(values, bins=60, color="#60a5fa", edgecolor="white")
        fixed_lower, fixed_upper = fixed_limits[column]
        if np.isfinite(fixed_lower):
            ax.axvline(fixed_lower, color="#059669", linestyle=":", label="Fixed minimum")
        if np.isfinite(fixed_upper):
            ax.axvline(fixed_upper, color="#059669", linestyle=":", label="Fixed maximum")
        if np.isfinite(lower):
            ax.axvline(lower, color="#dc2626", linestyle="--", label=f"Lower {mad_counts[column]:g} MAD")
        if np.isfinite(upper):
            ax.axvline(upper, color="#dc2626", linestyle="--", label=f"Upper {mad_counts[column]:g} MAD")
        ax.set_xlabel(label)
        ax.set_ylabel("Cells")
        if np.isfinite(lower) or np.isfinite(upper):
            ax.legend(frameon=False)
    fig.suptitle(f"{batch}: per-measurement-set RNA QC")
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def plot_scrublet_scores(scores, predicted, threshold, batch, outpath):
    """Plot per-measurement-set Scrublet scores and the resolved call threshold."""
    scores = np.asarray(scores, dtype=float)
    predicted = np.asarray(predicted, dtype=bool)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(scores[~predicted], bins=60, color="#60a5fa", alpha=0.85, label="Singlet calls")
    if predicted.any():
        ax.hist(scores[predicted], bins=60, color="#f97316", alpha=0.8, label="Doublet calls")
    if np.isfinite(threshold):
        ax.axvline(threshold, color="#7c3aed", linestyle="--", linewidth=1.5, label=f"Threshold = {threshold:.3g}")
    ax.set(xlabel="Scrublet doublet score", ylabel="Cells", title=f"{batch}: Scrublet doublet calls")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def run_scrublet_with_pca_fallback(
    counts, expected_doublet_rate, requested_n_prin_comps, adaptive_pca_fallback,
):
    """Run Scrublet, retrying only its explicit PCA-dimension failure."""
    import scrublet as scr

    def run(n_prin_comps):
        model = scr.Scrublet(
            counts,
            expected_doublet_rate=expected_doublet_rate,
            random_state=42,
        )
        scores, predicted = model.scrub_doublets(n_prin_comps=n_prin_comps)
        return model, scores, predicted

    try:
        model, scores, predicted = run(requested_n_prin_comps)
        return model, scores, predicted, requested_n_prin_comps, False
    except ValueError as error:
        match = re.search(
            r"n_components=\d+ must be between 1 and "
            r"min\(n_samples, n_features\)=(\d+)",
            str(error),
        )
        if not adaptive_pca_fallback or match is None:
            raise
        usable_dimension = int(match.group(1))
        fallback_n_prin_comps = min(requested_n_prin_comps - 1, usable_dimension - 1)
        if fallback_n_prin_comps < 1:
            raise ValueError(
                "Scrublet adaptive PCA fallback cannot run because fewer than "
                "two usable PCA dimensions remain"
            ) from error
        print(
            "Scrublet PCA fallback: requested "
            f"{requested_n_prin_comps}, retrying with {fallback_n_prin_comps} "
            f"because the usable dimension is {usable_dimension}."
        )
        model, scores, predicted = run(fallback_n_prin_comps)
        return model, scores, predicted, fallback_n_prin_comps, True


def prepare_matrix(adata, use_multimapping):
    if all(layer in adata.layers for layer in ("mature", "nascent", "ambiguous")):
        adata.X = (
            adata.layers["mature"].astype(np.float32)
            + adata.layers["nascent"].astype(np.float32)
            + adata.layers["ambiguous"].astype(np.float32)
        )
        for layer in ("mature", "nascent", "ambiguous"):
            adata.layers[layer] = normalize_sparse_index_dtypes(
                to_sparse_counts(adata.layers[layer])
            )
    if use_multimapping:
        if hasattr(adata.X, "data"):
            adata.X.data = np.round(adata.X.data)
        else:
            adata.X = np.round(adata.X)
    adata.X = normalize_sparse_index_dtypes(to_sparse_counts(adata.X))
    return adata


def main(args):
    if args.min_counts < 0:
        raise ValueError("Minimum RNA UMI counts must be non-negative")
    if not 0 < args.scrublet_expected_doublet_rate < 1:
        raise ValueError("Scrublet expected doublet rate must be in (0, 1)")
    mapping_dir = Path(args.mapping_dir)
    batch = batch_name(mapping_dir)
    label = safe_label(batch)
    folder = "counts_unfiltered_modified" if args.bc_replacement else "counts_unfiltered"
    adata = prepare_matrix(sc.read_h5ad(mapping_dir / folder / "adata.h5ad"), args.use_multimapping)

    symbols = pd.read_csv(mapping_dir / folder / "cells_x_genes.genes.names.txt", header=None)[0].astype(str)
    if len(symbols) != adata.n_vars:
        raise ValueError(f"{batch}: {len(symbols)} gene names for {adata.n_vars} variables")
    adata.var["symbol"] = symbols.to_numpy()
    adata.var_names = adata.var_names.astype(str).str.split(".").str[0]
    adata.var_names_make_unique()

    covariates = pd.read_csv(args.covariates, dtype=str)
    batch_column = covariates.columns[0]
    suffix = get_barcode_key(covariates, batch_column, batch)
    apply_batch_suffix(adata, batch, suffix, {}, record_corrected_barcode=args.bc_replacement)
    adata.obs[batch_column] = str(batch)
    covariates_for_obs = covariates.drop(columns=["barcode_key"], errors="ignore")
    adata.obs = adata.obs.join(covariates_for_obs.set_index(batch_column), on=batch_column, rsuffix="_covariate")
    adata.obs["batch_number"] = 1

    args.qc_dir.mkdir(parents=True, exist_ok=True)
    cell_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    knee_df = pd.DataFrame({"sum": cell_counts, "barcodes": adata.obs_names.astype(str)})
    knee_df = knee_df.sort_values("sum", ascending=False).reset_index(drop=True)
    knee_df["sum_log"] = np.log1p(knee_df["sum"])
    knee_df["rank"] = np.arange(1, len(knee_df) + 1)
    point_1, point_2 = get_elbow_knee_points(knee_df["rank"].to_numpy(), knee_df["sum_log"].to_numpy())

    selected_rank = None
    knee_threshold = np.nan
    keep_knee = np.ones(adata.n_obs, dtype=bool)
    if args.barcode_filter != "none":
        selected = point_1 if args.barcode_filter == "knee" else point_2
        if selected is not None:
            selected_rank = max(1, min(len(knee_df), int(round(selected[0]))))
            knee_threshold = float(knee_df.loc[selected_rank - 1, "sum"])
            keep_knee = cell_counts >= knee_threshold
    plot_barcode_rank(
        knee_df, point_1, point_2, selected_rank, batch,
        args.qc_dir / f"knee_plot_scRNA_{label}.png",
        barcode_filter=args.barcode_filter, selected_threshold=knee_threshold,
    )

    input_cells = adata.n_obs
    adata = adata[keep_knee].copy()
    post_knee_cells = adata.n_obs
    mt_prefix = "MT-" if args.reference == "human" else "Mt-"
    adata.var["mt"] = adata.var["symbol"].str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var["symbol"].str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo"],
        inplace=True,
        log1p=True,
        percent_top=tuple(value for value in (50, 100, 200, 500) if value <= adata.n_vars) or None,
    )

    limits = {
        "total_counts": mad_limits(adata.obs["log1p_total_counts"], args.mad_total_counts),
        "n_genes": mad_limits(adata.obs["log1p_n_genes_by_counts"], args.mad_n_genes),
    }
    flow_rows = []
    snapshots = []
    input_after_knee = np.ones(adata.n_obs, dtype=bool)
    current_keep = input_after_knee.copy()
    flow_rows.append({
        "measurement_set": batch,
        "step_order": 1,
        "filter_name": "barcode_filter",
        "filter_label": f"QC_barcode_filter = {args.barcode_filter}",
        "threshold": (
            f"rank ≤ {selected_rank:,}; RNA UMIs ≥ {knee_threshold:,.0f}"
            if selected_rank is not None else
            ("manual mode; no knee applied" if args.barcode_filter == "none" else "no valid knee; skipped")
        ),
        "applied": selected_rank is not None,
        "cells_before": input_cells,
        "cells_after": post_knee_cells,
        "cells_removed": input_cells - post_knee_cells,
        "removed_percent": 100 * (input_cells - post_knee_cells) / input_cells if input_cells else 0,
        "retained_percent_of_input": 100 * post_knee_cells / input_cells if input_cells else 0,
    })

    def apply_filter(
        filter_name, filter_label, threshold, values, lower, upper, enabled,
        bound_kind="fixed", lower_label=None, upper_label=None,
    ):
        nonlocal current_keep
        before_keep = current_keep.copy()
        candidate = np.ones(adata.n_obs, dtype=bool)
        if enabled:
            candidate &= values >= lower
            candidate &= values <= upper
        current_keep &= candidate
        before_count = int(before_keep.sum())
        after_count = int(current_keep.sum())
        flow_rows.append({
            "measurement_set": batch,
            "step_order": len(flow_rows) + 1,
            "filter_name": filter_name,
            "filter_label": filter_label,
            "threshold": threshold,
            "applied": bool(enabled),
            "cells_before": before_count,
            "cells_after": after_count,
            "cells_removed": before_count - after_count,
            "removed_percent": 100 * (before_count - after_count) / before_count if before_count else 0,
            "retained_percent_of_input": 100 * after_count / input_cells if input_cells else 0,
        })
        snapshots.append({
            "filter_label": filter_label,
            "metric_label": filter_name.replace("_", " "),
            "values_before": values[before_keep],
            "values_after": values[current_keep],
            "lower": lower,
            "upper": upper,
            "bound_kind": bound_kind,
            "lower_label": lower_label,
            "upper_label": upper_label,
            "applied": bool(enabled),
            "removed": before_count - after_count,
            "cells_before": before_count,
            "cells_after": after_count,
        })

    total_counts = adata.obs["total_counts"].to_numpy(dtype=float)
    n_genes = adata.obs["n_genes_by_counts"].to_numpy(dtype=float)
    apply_filter(
        "total RNA UMIs", "QC_min_counts_per_cell",
        f"total_counts ≥ {args.min_counts:,}", total_counts, float(args.min_counts), np.inf,
        args.min_counts > 0, lower_label=f"Fixed minimum ({args.min_counts:,})",
    )
    non_mito_mad_specs = (
        ("log1p total RNA UMIs", "QC_MAD_total_counts", "total_counts", "log1p_total_counts", args.mad_total_counts, False),
        ("log1p detected genes", "QC_MAD_n_genes", "n_genes", "log1p_n_genes_by_counts", args.mad_n_genes, False),
    )
    for metric_name, parameter, key, column, n_mads, upper_only in non_mito_mad_specs:
        _median, metric_mad, lower, upper = limits[key]
        effective = n_mads > 0 and np.isfinite(metric_mad) and metric_mad > 0
        bound_text = (
            f"≤ {upper:.3g} ({n_mads:g} MAD upper tail)" if upper_only else
            f"{lower:.3g} ≤ value ≤ {upper:.3g} ({n_mads:g} MAD)"
        ) if effective else f"disabled (configured {n_mads:g}; observed MAD {metric_mad:.3g})"
        apply_filter(
            metric_name, parameter, bound_text,
            adata.obs[column].to_numpy(dtype=float), lower, upper, effective,
            bound_kind="MAD", lower_label=f"Lower {n_mads:g} MAD", upper_label=f"Upper {n_mads:g} MAD",
        )

    pre_scrublet_keep = current_keep.copy()
    scrublet_removed = 0
    scrublet_status = "disabled"
    scrublet_skip_reason = ""
    scrublet_n_prin_comps_used = np.nan
    scrublet_pca_fallback_used = False
    if args.enable_scrublet:
        selected = np.flatnonzero(current_keep)
        try:
            scrub, scores, predicted, scrublet_n_prin_comps_used, scrublet_pca_fallback_used = run_scrublet_with_pca_fallback(
                adata.X[selected],
                args.scrublet_expected_doublet_rate,
                args.scrublet_n_prin_comps,
                args.scrublet_adaptive_pca_fallback,
            )
        except Exception as error:
            if args.scrublet_failure_policy == "error":
                raise
            scrublet_status = "skipped_error"
            scrublet_skip_reason = " ".join(
                f"{type(error).__name__}: {error}".split()
            )
            before_count = int(pre_scrublet_keep.sum())
            print(
                f"{batch}: Scrublet skipped by SCRUBLET_failure_policy=skip "
                f"after {scrublet_skip_reason}"
            )
            flow_rows.append({
                "measurement_set": batch,
                "step_order": len(flow_rows) + 1,
                "filter_name": "scrublet_doublet_score",
                "filter_label": "Scrublet doublet removal",
                "threshold": f"skipped after error: {scrublet_skip_reason}",
                "applied": False,
                "cells_before": before_count,
                "cells_after": before_count,
                "cells_removed": 0,
                "removed_percent": 0,
                "retained_percent_of_input": 100 * before_count / input_cells if input_cells else 0,
            })
        else:
            scrublet_status = "applied"
            threshold = float(scrub.threshold_) if scrub.threshold_ is not None else np.nan
            adata.obs["doublet_scores"] = np.nan
            adata.obs.iloc[selected, adata.obs.columns.get_loc("doublet_scores")] = scores
            adata.obs["predicted_doublets"] = False
            predicted_column = adata.obs.columns.get_loc("predicted_doublets")
            adata.obs.iloc[selected, predicted_column] = predicted
            adata.obs["doublet_info"] = adata.obs["predicted_doublets"].astype(str)
            current_keep[selected[predicted]] = False
            scrublet_removed = int(predicted.sum())
            before_count = int(pre_scrublet_keep.sum())
            after_count = int(current_keep.sum())
            flow_rows.append({
                "measurement_set": batch,
                "step_order": len(flow_rows) + 1,
                "filter_name": "scrublet_doublet_score",
                "filter_label": "Scrublet doublet removal",
                "threshold": (
                    f"expected_doublet_rate = {args.scrublet_expected_doublet_rate:g}; "
                    f"PCA components = {scrublet_n_prin_comps_used}"
                    f"{' (adaptive fallback)' if scrublet_pca_fallback_used else ''}; "
                    f"call threshold = {threshold:.3g}"
                ),
                "applied": True,
                "cells_before": before_count,
                "cells_after": after_count,
                "cells_removed": scrublet_removed,
                "removed_percent": 100 * scrublet_removed / before_count if before_count else 0,
                "retained_percent_of_input": 100 * after_count / input_cells if input_cells else 0,
            })
            snapshots.append({
                "filter_label": "Scrublet doublet removal",
                "metric_label": "Scrublet doublet score",
                "values_before": scores,
                "values_after": scores[~predicted],
                "lower": -np.inf,
                "upper": threshold,
                "bound_kind": "Scrublet",
                "lower_label": None,
                "upper_label": f"Scrublet threshold ({threshold:.3g})",
                "applied": True,
                "removed": scrublet_removed,
                "cells_before": before_count,
                "cells_after": after_count,
            })
            plot_scrublet_scores(
                scores, predicted, threshold, batch,
                args.qc_dir / f"scrublet_scores_scRNA_{label}.png",
            )
    else:
        count = int(current_keep.sum())
        flow_rows.append({
            "measurement_set": batch,
            "step_order": len(flow_rows) + 1,
            "filter_name": "scrublet_doublet_score",
            "filter_label": "Scrublet doublet removal",
            "threshold": "disabled",
            "applied": False,
            "cells_before": count,
            "cells_after": count,
            "cells_removed": 0,
            "removed_percent": 0,
            "retained_percent_of_input": 100 * count / input_cells if input_cells else 0,
        })
    keep = current_keep

    fixed_plot_limits = {
        "log1p_total_counts": (
            np.log1p(args.min_counts) if args.min_counts > 0 else -np.inf,
            np.inf,
        ),
        "log1p_n_genes_by_counts": (-np.inf, np.inf),
        "pct_counts_mt": (-np.inf, np.inf),
    }
    plot_qc_distributions(
        adata.obs,
        limits,
        fixed_plot_limits,
        {
            "log1p_total_counts": args.mad_total_counts,
            "log1p_n_genes_by_counts": args.mad_n_genes,
            "pct_counts_mt": 0,
        },
        batch,
        args.qc_dir / f"qc_distributions_scRNA_{label}.png",
    )
    flow = pd.DataFrame(flow_rows)
    flow.to_csv(args.qc_dir / f"rna_qc_filter_flow_{label}.tsv", sep="\t", index=False)
    plot_filter_flow(flow, batch, args.qc_dir / f"rna_qc_filter_flow_{label}.png")
    plot_filter_steps(
        snapshots,
        batch,
        args.qc_dir / f"rna_qc_filter_steps_{label}.png",
        {
            "knee_df": knee_df,
            "selected_rank": selected_rank,
            "selected_threshold": knee_threshold,
            "cells_before": input_cells,
            "cells_after": post_knee_cells,
            "barcode_filter": args.barcode_filter,
        },
    )
    retained = adata[keep].copy()
    retained.write_h5ad(f"{label}_filtered.h5ad")

    row = {
        "measurement_set": batch,
        "input_barcodes": input_cells,
        "post_knee_cells": post_knee_cells,
        "post_min_counts_cells": int((total_counts >= args.min_counts).sum()) if args.min_counts > 0 else post_knee_cells,
        "post_mad_cells": int(pre_scrublet_keep.sum()),
        "retained_cells": retained.n_obs,
        "retained_fraction": retained.n_obs / input_cells if input_cells else np.nan,
        "removed_by_scrublet": scrublet_removed,
        "barcode_filter": args.barcode_filter,
        "knee_rank": selected_rank,
        "knee_umi_threshold": knee_threshold,
        "fixed_min_counts": args.min_counts,
        "mad_total_counts_n": args.mad_total_counts,
        "mad_n_genes_n": args.mad_n_genes,
        "scrublet_enabled": args.enable_scrublet,
        "scrublet_failure_policy": args.scrublet_failure_policy,
        "scrublet_status": scrublet_status,
        "scrublet_skip_reason": scrublet_skip_reason,
        "scrublet_expected_doublet_rate": args.scrublet_expected_doublet_rate,
        "scrublet_n_prin_comps_requested": args.scrublet_n_prin_comps,
        "scrublet_n_prin_comps_used": scrublet_n_prin_comps_used,
        "scrublet_pca_fallback_enabled": args.scrublet_adaptive_pca_fallback,
        "scrublet_pca_fallback_used": scrublet_pca_fallback_used,
    }
    for key, prefix in (("total_counts", "total_counts"), ("n_genes", "n_genes")):
        median, mad, lower, upper = limits[key]
        row.update({
            f"{prefix}_median": median,
            f"{prefix}_mad": mad,
            f"{prefix}_mad_lower": lower if np.isfinite(lower) else np.nan,
            f"{prefix}_mad_upper": upper if np.isfinite(upper) else np.nan,
        })
    pd.DataFrame([row]).to_csv(args.qc_dir / f"measurement_set_qc_{label}.tsv", sep="\t", index=False)
    print(f"{batch}: retained {retained.n_obs}/{input_cells} cells after per-measurement-set QC")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mapping_dir")
    parser.add_argument("covariates")
    parser.add_argument("--qc-dir", type=Path, required=True)
    parser.add_argument("--min-counts", type=int, default=500)
    parser.add_argument("--reference", choices=["human", "mouse"], required=True)
    parser.add_argument("--barcode-filter", choices=["none", "knee", "knee2"], default="knee")
    parser.add_argument("--mad-total-counts", type=float, default=5)
    parser.add_argument("--mad-n-genes", type=float, default=5)
    parser.add_argument("--enable-scrublet", action="store_true")
    parser.add_argument("--scrublet-expected-doublet-rate", type=float, default=0.08)
    parser.add_argument("--scrublet-n-prin-comps", type=int, default=30)
    parser.add_argument("--scrublet-adaptive-pca-fallback", action="store_true")
    parser.add_argument("--scrublet-failure-policy", choices=["error", "skip"], default="error")
    parser.add_argument("--bc-replacement", action="store_true")
    parser.add_argument("--use-multimapping", action="store_true")
    parsed = parser.parse_args()
    for name in ("mad_total_counts", "mad_n_genes"):
        if getattr(parsed, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be non-negative")
    if parsed.scrublet_n_prin_comps < 1:
        parser.error("--scrublet-n-prin-comps must be at least 1")
    main(parsed)
