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


def plot_barcode_rank(knee_df, point_1, point_2, selected_rank, batch, outpath):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(knee_df["rank"], knee_df["sum_log"], linewidth=1.2, color="#2563eb")
    for point, color, label in ((point_1, "#dc2626", "Knee 1"), (point_2, "#f59e0b", "Knee 2")):
        if point is not None:
            ax.axvline(int(round(point[0])), color=color, linestyle="--", label=label)
    if selected_rank is not None:
        ax.axvline(selected_rank, color="#111827", linestyle=":", label="Selected")
    ax.set(xlabel="Barcode rank", ylabel="Log1p RNA UMI counts", title=f"{batch}: barcode-rank knee")
    if point_1 is not None or point_2 is not None or selected_rank is not None:
        ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outpath, dpi=180, facecolor="white")
    plt.close(fig)


def plot_qc_distributions(obs, limits, batch, outpath):
    specs = [
        ("log1p_total_counts", "Log1p total RNA UMIs", limits["total_counts"]),
        ("log1p_n_genes_by_counts", "Log1p detected genes", limits["n_genes"]),
        ("pct_counts_mt", "Mitochondrial counts (%)", limits["pct_mito"]),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, (column, label, (_median, _mad, lower, upper)) in zip(axes, specs):
        values = pd.to_numeric(obs[column], errors="coerce").dropna()
        ax.hist(values, bins=60, color="#60a5fa", edgecolor="white")
        if np.isfinite(lower):
            ax.axvline(lower, color="#dc2626", linestyle="--", label="MAD bound")
        if np.isfinite(upper):
            ax.axvline(upper, color="#dc2626", linestyle="--", label="MAD bound")
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
    plot_barcode_rank(knee_df, point_1, point_2, selected_rank, batch, args.qc_dir / f"knee_plot_scRNA_{label}.png")

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
        percent_top=(50, 100, 200, 500),
    )

    limits = {
        "total_counts": mad_limits(adata.obs["log1p_total_counts"], args.mad_total_counts),
        "n_genes": mad_limits(adata.obs["log1p_n_genes_by_counts"], args.mad_n_genes),
        "pct_mito": mad_limits(adata.obs["pct_counts_mt"], args.mad_pct_mito, upper_only=True),
    }
    fixed_keep = np.ones(adata.n_obs, dtype=bool)
    fixed_keep &= adata.obs["total_counts"].to_numpy() >= args.min_counts
    if args.barcode_filter == "none":
        fixed_keep &= adata.obs["n_genes_by_counts"].to_numpy() >= args.min_genes
    fixed_keep &= adata.obs["pct_counts_mt"].to_numpy() < args.pct_mito
    mad_keep = np.ones(adata.n_obs, dtype=bool)
    for key, column in (("total_counts", "log1p_total_counts"), ("n_genes", "log1p_n_genes_by_counts"), ("pct_mito", "pct_counts_mt")):
        _median, _mad, lower, upper = limits[key]
        values = adata.obs[column].to_numpy(dtype=float)
        mad_keep &= (values >= lower) & (values <= upper)
    keep = fixed_keep & mad_keep

    plot_qc_distributions(adata.obs, limits, batch, args.qc_dir / f"qc_distributions_scRNA_{label}.png")
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
