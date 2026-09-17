#!/usr/bin/env python3
"""Apply HTO support and singlet filters after guide/clone cell QC."""

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


NEGATIVE_LABELS = {"negative", "doublet", "multiplet", "multiplets", "nan", "none", ""}


def safe_label(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "unknown"


def normalized_hto_classes(obs):
    if "hto_type_split" in obs:
        labels = obs["hto_type_split"].astype(str).str.strip()
    elif "hto_type" in obs:
        raw = obs["hto_type"].astype(str).str.strip()
        labels = raw.map(lambda value: "multiplets" if "-" in value or "," in value else value)
    else:
        raise ValueError("Hashing modality lacks hto_type_split/hto_type demultiplex annotations")
    lowered = labels.str.lower()
    labels = labels.mask(lowered.eq("multiplet"), "multiplets")
    return labels


def plot_batch(frame, flow, minimum, label, outpath):
    positive = frame.loc[frame["positive_singlet"], "hto_label"].value_counts().sort_index()
    final = frame.loc[frame["retained"], "hto_label"].value_counts().reindex(positive.index, fill_value=0)
    names = positive.index.tolist()
    x = np.arange(len(names))
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    axes[0].bar(x, positive.to_numpy(), color="#60a5fa")
    axes[0].axhline(minimum, color="#dc2626", linestyle="--", label=f"Minimum positive cells = {minimum}")
    axes[0].set(title="Positive singlets before HTO support filter", ylabel="Cells", xlabel="HTO")
    axes[0].set_xticks(x, names, rotation=45, ha="right")
    axes[0].legend(frameon=False)
    axes[1].bar(x, final.to_numpy(), color="#34d399")
    axes[1].axhline(minimum, color="#dc2626", linestyle="--", label=f"Minimum positive cells = {minimum}")
    axes[1].set(title="Retained HTO singlets", ylabel="Cells", xlabel="HTO")
    axes[1].set_xticks(x, names, rotation=45, ha="right")
    axes[1].legend(frameon=False)
    axes[2].axis("off")
    y_positions = np.linspace(0.82, 0.18, len(flow))
    for index, (y, row) in enumerate(zip(y_positions, flow)):
        text = (
            f"Step {row['step_order']}: {row['filter_label']}\n{row['threshold']}\n"
            f"{row['cells_before']:,} → {row['cells_after']:,} cells "
            f"({row['cells_removed']:,} removed)"
        )
        axes[2].text(0.5, y, text, ha="center", va="center", fontsize=10,
                     bbox={"boxstyle": "round,pad=0.45", "facecolor": "#f8fafc", "edgecolor": "#94a3b8"})
        if index + 1 < len(flow):
            axes[2].annotate("", xy=(0.5, y_positions[index + 1] + 0.10), xytext=(0.5, y - 0.10),
                             arrowprops={"arrowstyle": "->", "color": "#64748b"})
    fig.suptitle(f"{label}: post-clone HTO filtering", fontsize=14)
    fig.tight_layout()
    fig.savefig(outpath, dpi=170, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    import mudata as md

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_mudata")
    parser.add_argument("output_mudata")
    parser.add_argument("--outdir", type=Path, default=Path("hto_qc"))
    parser.add_argument("--min-positive-cells", type=int, default=20)
    parser.add_argument("--singlet-only", choices=["true", "false"], default="true")
    parser.add_argument("--batch-column", default="batch")
    parser.add_argument("--filtered-hashing-output", default="post_clone_hashing_filtered.h5ad")
    parser.add_argument("--unfiltered-hashing-output", default="post_clone_hashing_unfiltered.h5ad")
    args = parser.parse_args()
    if args.min_positive_cells < 1:
        parser.error("--min-positive-cells must be >= 1")
    singlet_only = args.singlet_only == "true"

    args.outdir.mkdir(parents=True, exist_ok=True)
    mdata = md.read_h5mu(args.input_mudata)
    if "hashing" not in mdata.mod:
        raise ValueError("MuData has no hashing modality")
    hashing = mdata.mod["hashing"]
    if not hashing.obs_names.equals(mdata.obs_names):
        raise ValueError("Hashing and MuData cell orders differ; refusing unsafe HTO filtering")
    if args.batch_column in hashing.obs:
        batches = hashing.obs[args.batch_column].astype(str).fillna("unknown")
        resolved_batch_column = args.batch_column
    elif "batch" in hashing.obs:
        batches = hashing.obs["batch"].astype(str).fillna("unknown")
        resolved_batch_column = "batch"
    else:
        batches = pd.Series("all", index=hashing.obs_names, dtype=str)
        resolved_batch_column = "unavailable"

    labels = normalized_hto_classes(hashing.obs)
    positive_singlet = ~labels.str.lower().isin(NEGATIVE_LABELS)
    frame = pd.DataFrame({
        "cell_barcode": hashing.obs_names.astype(str),
        "measurement_set": batches.to_numpy(),
        "hto_label": labels.to_numpy(),
        "positive_singlet": positive_singlet.to_numpy(),
    })
    support = (
        frame.loc[frame["positive_singlet"]]
        .groupby(["measurement_set", "hto_label"], observed=True)
        .size().rename("positive_cells").reset_index()
    )
    frame = frame.merge(support, how="left", on=["measurement_set", "hto_label"])
    frame["positive_cells"] = frame["positive_cells"].fillna(0).astype(int)
    frame["hto_called"] = frame["positive_singlet"] & (frame["positive_cells"] >= args.min_positive_cells)
    keep_support = ~frame["positive_singlet"] | frame["hto_called"]
    keep_singlet = frame["positive_singlet"] if singlet_only else pd.Series(True, index=frame.index)
    frame["retained"] = keep_support & keep_singlet
    frame["filter_reason"] = np.select(
        [~keep_support, keep_support & ~keep_singlet],
        ["hto_below_minimum_positive_cells", "not_positive_singlet"],
        default="retained",
    )
    frame.to_csv(args.outdir / "hto_cell_filter.tsv", sep="\t", index=False)
    support["called"] = support["positive_cells"] >= args.min_positive_cells
    support.to_csv(args.outdir / "hto_positive_cell_support.tsv", sep="\t", index=False)

    flows = []
    for label in sorted(frame["measurement_set"].unique()):
        selected = frame["measurement_set"] == label
        before = int(selected.sum())
        after_support = int((selected & keep_support).sum())
        after_singlet = int((selected & frame["retained"]).sum())
        batch_flow = [
            {
                "measurement_set": label, "step_order": 1,
                "filter_label": "HTO_min_positive_cells",
                "threshold": f"HTO called when positive singlet cells >= {args.min_positive_cells}",
                "applied": True, "cells_before": before, "cells_after": after_support,
                "cells_removed": before - after_support,
            },
            {
                "measurement_set": label, "step_order": 2,
                "filter_label": "HTO_keep_singlets_only",
                "threshold": "retain positive HTO singlets only" if singlet_only else "disabled",
                "applied": singlet_only, "cells_before": after_support, "cells_after": after_singlet,
                "cells_removed": after_support - after_singlet,
            },
        ]
        for row in batch_flow:
            row["removed_percent"] = 100 * row["cells_removed"] / row["cells_before"] if row["cells_before"] else 0.0
            row["retained_percent_of_input"] = 100 * row["cells_after"] / before if before else 0.0
            row["batch_column"] = resolved_batch_column
        flows.extend(batch_flow)
        plot_batch(
            frame.loc[selected], batch_flow, args.min_positive_cells, label,
            args.outdir / f"hto_filter_steps_{safe_label(label)}.png",
        )
    pd.DataFrame(flows).to_csv(args.outdir / "hto_filter_flow.tsv", sep="\t", index=False)

    keep = frame["retained"].to_numpy(dtype=bool)
    metrics = {
        "execution_order": "after_guide_assignment_and_clone_removal",
        "min_positive_cells": args.min_positive_cells,
        "singlet_only": singlet_only,
        "batch_column": resolved_batch_column,
        "input_cells": int(len(frame)),
        "retained_cells": int(keep.sum()),
        "removed_cells": int((~keep).sum()),
        "called_htos": int(support["called"].sum()),
    }
    (args.outdir / "hto_filter_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    hashing.obs["hto_positive_cells_post_clone"] = frame["positive_cells"].to_numpy(dtype=int)
    hashing.obs["hto_called_post_clone"] = frame["hto_called"].to_numpy(dtype=bool)
    hashing.write_h5ad(args.unfiltered_hashing_output)
    hashing[keep].copy().write_h5ad(args.filtered_hashing_output)
    filtered = mdata[keep].copy()
    filtered.uns["hto_post_clone_filter"] = metrics
    filtered.write_h5mu(args.output_mudata)
    print(f"Post-clone HTO QC retained {int(keep.sum())}/{len(keep)} cells")


if __name__ == "__main__":
    main()
