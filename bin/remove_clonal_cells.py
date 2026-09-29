#!/usr/bin/env python3
"""Detect guide-barcode clones and optionally remove their cells.

The grouping test follows Wang et al. (BMC Genomics 2022) and the MIT-licensed
reference implementation at https://github.com/yihan1119/Group_clone.  Cells
are visited deterministically and their assigned-guide overlap is compared
with the representative cell of each previously discovered clone using a
Bonferroni-corrected hypergeometric survival probability.  A cell matching
more than one clone is marked as an ambiguous clone/doublet.
"""

import argparse
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import hypergeom


def _assigned_matrix(guide):
    if "guide_assignment" not in guide.layers:
        raise ValueError("guide.layers['guide_assignment'] is required for clone detection")
    matrix = guide.layers["guide_assignment"]
    matrix = matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)
    matrix = matrix.astype(bool).astype(np.uint8)
    matrix.eliminate_zeros()
    return matrix


def group_clones(assignment, alpha=0.05):
    """Faithful, deterministic implementation of the published algorithm."""
    n_cells, n_guides = assignment.shape
    threshold = alpha / max(1.0, n_cells * n_cells / 2.0)
    guide_counts = np.asarray(assignment.sum(axis=1)).ravel().astype(int)
    representatives = []
    members = []
    ambiguous = np.zeros(n_cells, dtype=bool)
    clone_index = np.full(n_cells, -1, dtype=int)
    representatives_by_guide = defaultdict(list)

    for cell_idx in range(n_cells):
        guides = assignment.indices[assignment.indptr[cell_idx] : assignment.indptr[cell_idx + 1]]
        candidate_ids = sorted(
            {clone_id for guide_idx in guides for clone_id in representatives_by_guide[int(guide_idx)]}
        )
        matches = []
        if candidate_ids:
            rep_rows = np.asarray([representatives[c] for c in candidate_ids], dtype=int)
            overlaps = np.asarray(assignment[rep_rows].dot(assignment[cell_idx].T).todense()).ravel()
            rep_guide_counts = guide_counts[rep_rows]
            pvalues = hypergeom.sf(
                overlaps - 1,
                n_guides,
                rep_guide_counts,
                guide_counts[cell_idx],
            )
            matches = [candidate_ids[i] for i in np.flatnonzero(pvalues <= threshold)]

        if len(matches) == 0:
            clone_id = len(representatives)
            representatives.append(cell_idx)
            members.append([cell_idx])
            clone_index[cell_idx] = clone_id
            for guide_idx in guides:
                representatives_by_guide[int(guide_idx)].append(clone_id)
        elif len(matches) == 1:
            clone_id = matches[0]
            members[clone_id].append(cell_idx)
            clone_index[cell_idx] = clone_id
        else:
            ambiguous[cell_idx] = True

    return {
        "members": members,
        "representatives": np.asarray(representatives, dtype=int),
        "clone_index": clone_index,
        "ambiguous": ambiguous,
        "guide_counts": guide_counts,
        "bonferroni_threshold": threshold,
    }


def _gene_totals(mdata):
    gene = mdata.mod["gene"]
    if "total_gene_umis" in gene.obs:
        return gene.obs["total_gene_umis"].to_numpy(dtype=float)
    return np.asarray(gene.X.sum(axis=1)).ravel().astype(float)


def decide_cells(grouped, action, min_clone_size, gene_totals):
    n_cells = len(grouped["clone_index"])
    clone_sizes = np.ones(n_cells, dtype=int)
    representative = np.zeros(n_cells, dtype=bool)
    clonal = np.zeros(n_cells, dtype=bool)
    for clone_id, member_list in enumerate(grouped["members"]):
        size = len(member_list)
        clone_sizes[member_list] = size
        if size >= min_clone_size:
            clonal[member_list] = True
        best = member_list[int(np.argmax(gene_totals[member_list]))]
        representative[best] = True

    keep = np.ones(n_cells, dtype=bool)
    keep[grouped["ambiguous"]] = False
    if action == "drop_clonal":
        keep[clonal] = False
    elif action == "keep_representative":
        keep[clonal & ~representative] = False
    elif action != "mark_only":
        raise ValueError(f"Unknown clone-removal action: {action}")
    if action == "mark_only":
        keep[:] = True
    return keep, clonal, clone_sizes, representative


def _plots(assignments, clone_summary, outdir, applicability):
    sizes = clone_summary.loc[clone_summary["clone_size"] >= 2, "clone_size"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    if sizes.size:
        bins = np.arange(1.5, sizes.max() + 1.6)
        axes[0].hist(sizes, bins=bins, color="#2878b5", edgecolor="white")
        axes[0].set_yscale("log")
    else:
        axes[0].text(0.5, 0.5, "No multi-cell clones", ha="center", va="center")
    axes[0].set_xlabel("Cells per detected clone")
    axes[0].set_ylabel("Clone count (log scale)")
    axes[0].set_title("Clone-size distribution")

    counts = assignments["status"].value_counts().reindex(
        ["non_clonal", "clonal_retained", "clonal_removed", "ambiguous_retained", "ambiguous_removed"], fill_value=0
    )
    axes[1].bar(counts.index, counts.values, color=["#6baed6", "#74c476", "#fb6a4a", "#bcbddc", "#756bb1"])
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].set_ylabel("Cells")
    axes[1].set_title("Clone-filter outcome")
    if applicability == "low_power_warning":
        fig.suptitle(
            "Low-power design for guide-barcode clone detection (<1,000 guides or median <10 guides/cell)",
            color="#b45309",
            fontsize=11,
        )
    fig.tight_layout(rect=(0, 0, 1, 0.94) if applicability == "low_power_warning" else None)
    fig.savefig(outdir / "clone_filter_summary.png", dpi=180, facecolor="white")
    plt.close(fig)


def main():
    import mudata as md

    parser = argparse.ArgumentParser()
    parser.add_argument("input_mudata")
    parser.add_argument("output_mudata")
    parser.add_argument("--outdir", default="clone_qc")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-clone-size", type=int, default=2)
    parser.add_argument(
        "--action",
        choices=["mark_only", "drop_clonal", "keep_representative"],
        default="drop_clonal",
    )
    args = parser.parse_args()
    if not 0 < args.alpha <= 1:
        parser.error("--alpha must be in (0, 1]")
    if args.min_clone_size < 2:
        parser.error("--min-clone-size must be at least 2")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mdata = md.read_h5mu(args.input_mudata)
    if "guide" not in mdata.mod or "gene" not in mdata.mod:
        raise ValueError("MuData must contain gene and guide modalities")
    guide = mdata.mod["guide"]
    if not guide.obs_names.equals(mdata.obs_names):
        raise ValueError("Guide and MuData cell orders differ; refusing unsafe clone filtering")

    assignment = _assigned_matrix(guide)
    grouped = group_clones(assignment, args.alpha)
    keep, clonal, clone_sizes, representative = decide_cells(
        grouped, args.action, args.min_clone_size, _gene_totals(mdata)
    )
    clone_labels = np.where(
        grouped["clone_index"] >= 0,
        np.char.add("clone_", grouped["clone_index"].astype(str)),
        "ambiguous",
    )
    status = np.full(mdata.n_obs, "non_clonal", dtype=object)
    status[clonal & keep] = "clonal_retained"
    status[clonal & ~keep] = "clonal_removed"
    status[grouped["ambiguous"] & keep] = "ambiguous_retained"
    status[grouped["ambiguous"] & ~keep] = "ambiguous_removed"
    audit = pd.DataFrame(
        {
            "cell_barcode": mdata.obs_names.astype(str),
            "clone_id": clone_labels,
            "clone_size": clone_sizes,
            "assigned_guides": grouped["guide_counts"],
            "clone_representative": representative,
            "ambiguous_clone_match": grouped["ambiguous"],
            "retained": keep,
            "status": status,
        }
    )
    audit.to_csv(outdir / "clone_cell_assignments.tsv", sep="\t", index=False)
    clone_summary = (
        audit.loc[~audit["ambiguous_clone_match"]]
        .groupby("clone_id", observed=True)
        .agg(clone_size=("cell_barcode", "size"), retained_cells=("retained", "sum"))
        .reset_index()
        .sort_values(["clone_size", "clone_id"], ascending=[False, True])
    )
    clone_summary.to_csv(outdir / "clone_groups.tsv", sep="\t", index=False)

    median_guides = float(np.median(grouped["guide_counts"]))
    applicability = "supported" if assignment.shape[1] >= 1000 and median_guides >= 10 else "low_power_warning"
    metrics = {
        "action": args.action,
        "alpha": args.alpha,
        "bonferroni_threshold": grouped["bonferroni_threshold"],
        "min_clone_size": args.min_clone_size,
        "input_cells": int(mdata.n_obs),
        "retained_cells": int(keep.sum()),
        "removed_cells": int((~keep).sum()),
        "guide_library_size": int(assignment.shape[1]),
        "median_assigned_guides_per_cell": median_guides,
        "mean_assigned_guides_per_cell": float(np.mean(grouped["guide_counts"])),
        "detected_clone_groups": int((clone_summary["clone_size"] >= args.min_clone_size).sum()),
        "cells_in_detected_clones": int(clonal.sum()),
        "ambiguous_clone_matches": int(grouped["ambiguous"].sum()),
        "clonal_cell_fraction": float(clonal.mean()),
        "applicability": applicability,
        "applicability_note": (
            "Wang et al. recommend at least 1,000 guides and at least 10 detected guides per cell."
        ),
        "method_source": "Wang et al., BMC Genomics 2022; github.com/yihan1119/Group_clone",
    }
    pd.DataFrame([metrics]).to_csv(outdir / "clone_metrics.tsv", sep="\t", index=False)
    (outdir / "clone_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    _plots(audit, clone_summary, outdir, applicability)

    filtered = mdata[keep].copy()
    retained_audit = audit.loc[keep].set_index("cell_barcode")
    for frame in [filtered.obs, filtered.mod["gene"].obs, filtered.mod["guide"].obs]:
        aligned = retained_audit.reindex(frame.index.astype(str))
        frame["clone_id"] = aligned["clone_id"].to_numpy()
        frame["clone_size"] = aligned["clone_size"].to_numpy(dtype=int)
        frame["clone_representative"] = aligned["clone_representative"].to_numpy(dtype=bool)
    filtered.uns["clone_removal"] = metrics
    filtered.write_h5mu(args.output_mudata)


if __name__ == "__main__":
    main()
