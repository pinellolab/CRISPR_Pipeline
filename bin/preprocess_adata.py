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
    for point, color, label in ((point_1, "#dc2626", "Knee 1"), (point_2, "#f59e0b", "Knee 2")):
        if point is not None:
            ax.axvline(int(round(point[0])), color=color, linestyle="--", label=label)
    if selected_rank is not None:
        ax.axvline(selected_rank, color="#111827", linestyle=":", label="Selected")
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
    ax.set(xlabel="Barcode rank", ylabel="Log1p RNA UMI counts", title=f"{batch}: barcode-rank knee")
    if point_1 is not None or point_2 is not None or selected_rank is not None:
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def plot_filter_steps(snapshots, batch, outpath):
    """Plot each cell filter immediately before and after it is applied."""
    fig, axes = plt.subplots(len(snapshots), 3, figsize=(17, 3.15 * len(snapshots)), squeeze=False)
    for row_index, snapshot in enumerate(snapshots):
        before = np.asarray(snapshot["values_before"], dtype=float)
        after = np.asarray(snapshot["values_after"], dtype=float)
        before = before[np.isfinite(before)]
        after = after[np.isfinite(after)]
        finite = before if before.size else np.array([0.0])
        low, high = np.nanpercentile(finite, [0.5, 99.5]) if finite.size > 1 else (finite[0] - 0.5, finite[0] + 0.5)
        if not np.isfinite(low) or not np.isfinite(high) or low == high:
            low, high = float(np.nanmin(finite)) - 0.5, float(np.nanmax(finite)) + 0.5
        bins = np.linspace(low, high, 51)
        for column, (values, color, state) in enumerate(
            ((before, "#60a5fa", "Before"), (after, "#34d399", "After"))
        ):
            ax = axes[row_index, column]
            ax.hist(values, bins=bins, color=color, edgecolor="white")
            add_filter_bounds(ax, snapshot)
            count_key = "cells_before" if state == "Before" else "cells_after"
            ax.set_title(f"{state}: {snapshot[count_key]:,} cells")
            ax.set_xlabel(snapshot["metric_label"])
            ax.set_ylabel("Cells")
        box_ax = axes[row_index, 2]
        box_values = [before if before.size else np.array([np.nan]), after if after.size else np.array([np.nan])]
        boxes = box_ax.boxplot(
            box_values, vert=False, labels=["Before", "After"], showfliers=False,
            patch_artist=True, medianprops={"color": "#111827", "linewidth": 1.5},
        )
        for patch, color in zip(boxes["boxes"], ("#93c5fd", "#6ee7b7")):
            patch.set_facecolor(color)
        add_filter_bounds(box_ax, snapshot, show_legend=True)
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


def add_filter_bounds(ax, snapshot, show_legend=False):
    """Draw explicitly named fixed or MAD limits on a QC axis."""
    if not snapshot["applied"]:
        return
    color = "#dc2626" if snapshot["bound_kind"] == "MAD" else "#059669"
    linestyle = "--" if snapshot["bound_kind"] == "MAD" else ":"
    for bound, label in (
        (snapshot["lower"], snapshot.get("lower_label")),
        (snapshot["upper"], snapshot.get("upper_label")),
    ):
        if np.isfinite(bound):
            ax.axvline(bound, color=color, linestyle=linestyle, linewidth=1.5, label=label)
    if show_legend and (snapshot.get("lower_label") or snapshot.get("upper_label")):
        ax.legend(frameon=False, fontsize=8, loc="best")


def plot_filter_flow(flow, batch, outpath):
    """Render cells → filter → cells in the exact execution order."""
    node_count = 1 + 2 * len(flow)
    fig, ax = plt.subplots(figsize=(11, max(10, node_count * 0.82)))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, node_count + 1)
    ax.axis("off")
    y = node_count
    initial = int(flow.iloc[0]["cells_before"]) if len(flow) else 0
    ax.text(0.5, y, f"{initial:,} input barcodes", ha="center", va="center", fontsize=12, fontweight="bold",
            bbox={"boxstyle": "round,pad=0.5", "facecolor": "#dbeafe", "edgecolor": "#2563eb"})
    for _, step in flow.iterrows():
        ax.annotate("", xy=(0.5, y - 0.72), xytext=(0.5, y - 0.28), arrowprops={"arrowstyle": "->", "color": "#64748b"})
        y -= 1
        status = "APPLIED" if bool(step["applied"]) else "DISABLED / SKIPPED"
        ax.text(
            0.5, y,
            f"{step['filter_label']}\n{step['threshold']}\n{status}",
            ha="center", va="center", fontsize=9.5,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f8fafc", "edgecolor": "#94a3b8"},
        )
        ax.annotate("", xy=(0.5, y - 0.72), xytext=(0.5, y - 0.28), arrowprops={"arrowstyle": "->", "color": "#64748b"})
        y -= 1
        removed = int(step["cells_removed"])
        after = int(step["cells_after"])
        color = "#dcfce7" if removed == 0 else "#fef3c7"
        ax.text(
            0.5, y,
            f"{after:,} cells retained  |  {removed:,} removed ({float(step['removed_percent']):.2f}%)",
            ha="center", va="center", fontsize=10, fontweight="bold",
            bbox={"boxstyle": "round,pad=0.45", "facecolor": color, "edgecolor": "#16a34a"},
        )
    fig.suptitle(f"{batch}: RNA QC filtering flow (pipeline order)", fontsize=15)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def plot_qc_distributions(obs, limits, fixed_limits, mad_counts, batch, outpath):
    specs = [
        ("log1p_total_counts", "Log1p total RNA UMIs", limits["total_counts"]),
        ("log1p_n_genes_by_counts", "Log1p detected genes", limits["n_genes"]),
        ("pct_counts_mt", "Mitochondrial counts (%)", limits["pct_mito"]),
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
    if args.min_genes < 0 or args.min_counts < 0:
        raise ValueError("Minimum genes and RNA UMI counts must be non-negative")
    if not 0 <= args.pct_mito <= 100:
        raise ValueError("Mitochondrial percentage cutoff must be in [0, 100]")
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
        "pct_mito": mad_limits(adata.obs["pct_counts_mt"], args.mad_pct_mito, upper_only=True),
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
    pct_mito = adata.obs["pct_counts_mt"].to_numpy(dtype=float)
    apply_filter(
        "total RNA UMIs", "QC_min_counts_per_cell",
        f"total_counts ≥ {args.min_counts:,}", total_counts, float(args.min_counts), np.inf,
        args.min_counts > 0, lower_label=f"Fixed minimum ({args.min_counts:,})",
    )
    apply_filter(
        "detected genes", "QC_min_genes_per_cell",
        (
            f"n_genes_by_counts ≥ {args.min_genes:,}"
            if args.barcode_filter == "none" else
            f"skipped because QC_barcode_filter = {args.barcode_filter}"
        ),
        n_genes, float(args.min_genes), np.inf,
        args.barcode_filter == "none" and args.min_genes > 0,
        lower_label=f"Fixed minimum ({args.min_genes:,})",
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

    # Mitochondrial filtering is deliberately last: the fixed cutoff followed by
    # the optional upper-tail MAD cutoff. This keeps its impact explicit in the
    # sequential audit and avoids conflating it with RNA-complexity filtering.
    apply_filter(
        "mitochondrial percentage", "QC_pct_mito",
        f"pct_counts_mt < {args.pct_mito:g}%", pct_mito, -np.inf,
        np.nextafter(float(args.pct_mito), -np.inf),
        args.pct_mito < 100, upper_label=f"Fixed maximum ({args.pct_mito:g}%)",
    )
    _median, mito_mad, lower, upper = limits["pct_mito"]
    mito_mad_effective = args.mad_pct_mito > 0 and np.isfinite(mito_mad) and mito_mad > 0
    mito_bound_text = (
        f"≤ {upper:.3g} ({args.mad_pct_mito:g} MAD upper tail)"
        if mito_mad_effective else
        f"disabled (configured {args.mad_pct_mito:g}; observed MAD {mito_mad:.3g})"
    )
    apply_filter(
        "mitochondrial percentage", "QC_MAD_pct_mito", mito_bound_text,
        pct_mito, lower, upper, mito_mad_effective,
        bound_kind="MAD", upper_label=f"Upper {args.mad_pct_mito:g} MAD",
    )
    keep = current_keep

    # Preserve the historical summary meanings even though the visual audit now
    # applies mitochondrial filters last. These masks are order-independent.
    fixed_keep = np.ones(adata.n_obs, dtype=bool)
    if args.min_counts > 0:
        fixed_keep &= total_counts >= args.min_counts
    if args.barcode_filter == "none" and args.min_genes > 0:
        fixed_keep &= n_genes >= args.min_genes
    if args.pct_mito < 100:
        fixed_keep &= pct_mito < args.pct_mito

    fixed_plot_limits = {
        "log1p_total_counts": (
            np.log1p(args.min_counts) if args.min_counts > 0 else -np.inf,
            np.inf,
        ),
        "log1p_n_genes_by_counts": (
            np.log1p(args.min_genes)
            if args.barcode_filter == "none" and args.min_genes > 0
            else -np.inf,
            np.inf,
        ),
        "pct_counts_mt": (-np.inf, args.pct_mito),
    }
    plot_qc_distributions(
        adata.obs,
        limits,
        fixed_plot_limits,
        {
            "log1p_total_counts": args.mad_total_counts,
            "log1p_n_genes_by_counts": args.mad_n_genes,
            "pct_counts_mt": args.mad_pct_mito,
        },
        batch,
        args.qc_dir / f"qc_distributions_scRNA_{label}.png",
    )
    flow = pd.DataFrame(flow_rows)
    flow.to_csv(args.qc_dir / f"rna_qc_filter_flow_{label}.tsv", sep="\t", index=False)
    plot_filter_flow(flow, batch, args.qc_dir / f"rna_qc_filter_flow_{label}.png")
    plot_filter_steps(snapshots, batch, args.qc_dir / f"rna_qc_filter_steps_{label}.png")
    retained = adata[keep].copy()
    retained.write_h5ad(f"{label}_filtered.h5ad")

    row = {
        "measurement_set": batch,
        "input_barcodes": input_cells,
        "post_knee_cells": post_knee_cells,
        "post_fixed_threshold_cells": int(fixed_keep.sum()),
        "retained_cells": retained.n_obs,
        "retained_fraction": retained.n_obs / input_cells if input_cells else np.nan,
        "removed_by_fixed_thresholds": int((~fixed_keep).sum()),
        "removed_by_mad_after_fixed": int(fixed_keep.sum() - keep.sum()),
        "barcode_filter": args.barcode_filter,
        "knee_rank": selected_rank,
        "knee_umi_threshold": knee_threshold,
        "fixed_min_counts": args.min_counts,
        "fixed_min_genes": args.min_genes if args.barcode_filter == "none" else np.nan,
        "fixed_pct_mito_max": args.pct_mito,
        "mad_total_counts_n": args.mad_total_counts,
        "mad_n_genes_n": args.mad_n_genes,
        "mad_pct_mito_n": args.mad_pct_mito,
    }
    for key, prefix in (("total_counts", "total_counts"), ("n_genes", "n_genes"), ("pct_mito", "pct_mito")):
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
    parser.add_argument("--min-genes", type=int, default=100)
    parser.add_argument("--min-counts", type=int, default=0)
    parser.add_argument("--pct-mito", type=float, default=20)
    parser.add_argument("--reference", choices=["human", "mouse"], required=True)
    parser.add_argument("--barcode-filter", choices=["none", "knee", "knee2"], default="knee")
    parser.add_argument("--mad-total-counts", type=float, default=0)
    parser.add_argument("--mad-n-genes", type=float, default=0)
    parser.add_argument("--mad-pct-mito", type=float, default=0)
    parser.add_argument("--bc-replacement", action="store_true")
    parser.add_argument("--use-multimapping", action="store_true")
    parsed = parser.parse_args()
    for name in ("mad_total_counts", "mad_n_genes", "mad_pct_mito"):
        if getattr(parsed, name) < 0:
            parser.error(f"--{name.replace('_', '-')} must be non-negative")
    main(parsed)
