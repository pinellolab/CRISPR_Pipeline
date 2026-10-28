#!/usr/bin/env python3
"""Filter cells with excessive assigned guides and create per-batch QC."""

import argparse
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp


def assigned_counts(guide):
    if "guide_assignment" not in guide.layers:
        raise ValueError("guide.layers['guide_assignment'] is required")
    matrix = guide.layers["guide_assignment"]
    if sp.issparse(matrix):
        return np.asarray((matrix > 0).sum(axis=1)).ravel().astype(int)
    return np.count_nonzero(np.asarray(matrix), axis=1).astype(int)


def safe(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "unknown"


def batch_values(mdata, requested):
    # MuData can expose common modality columns only as gene:batch/guide:batch.
    # Read aligned modality observations as well, preserving measurement-set QC.
    candidates = (requested, "batch", "measurement_sets", "measurement_set")
    for observations in (mdata.obs, mdata.mod["guide"].obs, mdata.mod["gene"].obs):
        if not observations.index.equals(mdata.obs_names):
            continue
        for column in candidates:
            if column and column in observations:
                return observations[column].astype("string").fillna("unknown").astype(str), column
    return pd.Series("all", index=mdata.obs_names, dtype=str), "unavailable"


def plot_filter(counts, keep, label, maximum, output, minimum=0):
    upper = max(int(counts.max()) if counts.size else 0, maximum, 1)
    bins = np.arange(-0.5, upper + 1.5, 1)
    before, after = len(counts), int(keep.sum())
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for ax, values, color, title in (
        (axes[0], counts, "#60a5fa", f"Before: {before:,} cells"),
        (axes[1], counts[keep], "#34d399", f"After: {after:,} cells"),
    ):
        ax.hist(values, bins=bins, color=color, edgecolor="white")
        if maximum > 0:
            ax.axvline(maximum + 0.5, color="#dc2626", linestyle="--", label=f"Maximum = {maximum}")
            ax.legend(frameon=False)
        ax.set(title=title, xlabel="Assigned gRNAs per cell", ylabel="Cells")
    axes[2].axis("off")
    text = (
        f"{before:,} cells before\n\n↓\n\n"
        + (f"{minimum} ≤ assigned gRNAs ≤ {maximum}\nAPPLIED" if maximum > 0 else f"assigned gRNAs ≥ {minimum}")
        + f"\n\n↓\n\n{after:,} cells after\n{before-after:,} removed"
    )
    axes[2].text(0.5, 0.5, text, ha="center", va="center", fontsize=11,
                 bbox={"boxstyle": "round,pad=0.6", "facecolor": "#f8fafc", "edgecolor": "#94a3b8"})
    fig.suptitle(f"{label}: guide-assignment cell filter", fontsize=14)
    fig.tight_layout()
    fig.savefig(output, dpi=170, facecolor="white")
    plt.close(fig)


def main():
    import mudata as md

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_mudata")
    parser.add_argument("output_mudata")
    parser.add_argument("--outdir", type=Path, default=Path("guide_assignment_qc"))
    parser.add_argument("--max-guides-per-cell", type=int, default=15)
    parser.add_argument("--min-guides-per-cell", type=int, default=1)
    parser.add_argument("--batch-column", default="batch")
    args = parser.parse_args()
    if args.max_guides_per_cell < 0:
        parser.error("--max-guides-per-cell must be >= 0; use 0 to disable")
    if args.min_guides_per_cell < 0 or (args.max_guides_per_cell and args.min_guides_per_cell > args.max_guides_per_cell):
        parser.error("Guide minimum must be nonnegative and cannot exceed a nonzero maximum")
    args.outdir.mkdir(parents=True, exist_ok=True)

    mdata = md.read_h5mu(args.input_mudata)
    guide = mdata.mod["guide"]
    if not guide.obs_names.equals(mdata.obs_names):
        raise ValueError("Guide and MuData cell orders differ")
    counts = assigned_counts(guide)
    keep = counts <= args.max_guides_per_cell if args.max_guides_per_cell else np.ones(len(counts), bool)
    keep &= counts >= args.min_guides_per_cell
    batches, batch_column = batch_values(mdata, args.batch_column)
    batch_array = batches.to_numpy()
    pd.DataFrame({
        "cell_barcode": mdata.obs_names.astype(str), "measurement_set": batch_array,
        "assigned_guides": counts, "retained": keep,
        "filter_reason": np.where(keep, "retained", np.where(counts < args.min_guides_per_cell,
                                "assigned_guides_below_minimum", "assigned_guides_above_maximum")),
    }).to_csv(args.outdir / "guide_assignment_cell_filter.tsv", sep="\t", index=False)

    rows = []
    groups = [("all", np.ones(len(counts), bool))]
    if batch_column != "unavailable":
        groups += [(label, batch_array == label) for label in sorted(np.unique(batch_array))]
    for label, selected in groups:
        before, after = int(selected.sum()), int(keep[selected].sum())
        removed = before - after
        rows.append({
            "measurement_set": label, "step_order": 1,
            "filter_label": "GUIDE_ASSIGNMENT_min_max_guides_per_cell",
            "threshold": f"minimum={args.min_guides_per_cell}; maximum={args.max_guides_per_cell or 'unlimited'}",
            "applied": args.max_guides_per_cell > 0 or args.min_guides_per_cell > 0, "cells_before": before,
            "cells_after": after, "cells_removed": removed,
            "removed_percent": 100 * removed / before if before else 0.0,
            "retained_percent_of_input": 100 * after / before if before else 0.0,
            "batch_column": batch_column,
        })
        plot_filter(counts[selected], keep[selected], label, args.max_guides_per_cell,
                    args.outdir / f"guide_assignment_filter_steps_{safe(label)}.png", args.min_guides_per_cell)
    pd.DataFrame(rows).to_csv(args.outdir / "guide_assignment_filter_flow.tsv", sep="\t", index=False)
    metrics = {
        "max_guides_per_cell": args.max_guides_per_cell, "min_guides_per_cell": args.min_guides_per_cell,
        "applied": args.max_guides_per_cell > 0 or args.min_guides_per_cell > 0,
        "batch_column": batch_column, "input_cells": len(keep),
        "retained_cells": int(keep.sum()), "removed_cells": int((~keep).sum()),
    }
    (args.outdir / "guide_assignment_filter_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    mdata.obs["assigned_guides_before_qc"] = counts
    for modality in mdata.mod.values():
        if modality.n_obs == mdata.n_obs and modality.obs_names.equals(mdata.obs_names):
            modality.obs["assigned_guides_before_qc"] = counts
    filtered = mdata[keep].copy()
    filtered.uns["guide_assignment_cell_filter"] = metrics
    filtered.write_h5mu(args.output_mudata)
    print(f"Guide-assignment QC retained {int(keep.sum())}/{len(keep)} cells")


if __name__ == "__main__":
    main()
