#!/usr/bin/env python3
"""Audit guide mapping, barcode overlap, and configured guide orientation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import mudata as md
import numpy as np
from scipy import sparse


def parse_bool(value: str) -> bool:
    value = str(value).strip().lower()
    if value in {"true", "1", "yes", "y"}:
        return True
    if value in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"not a boolean: {value}")


def read_expected_guides(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters="\t,")
        except csv.Error:
            dialect = csv.excel_tab
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"guide metadata has no header: {path}")
        id_column = next(
            (name for name in ("guide_id", "sgrna_id", "id", "name") if name in reader.fieldnames),
            None,
        )
        if id_column is None:
            raise ValueError(
                "guide metadata must contain one of guide_id, sgrna_id, id, or name; "
                f"found {reader.fieldnames}"
            )
        return [str(row[id_column]).strip() for row in reader if str(row.get(id_column, "")).strip()]


def counts_by_feature(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.sum(axis=0)).ravel()
    return np.asarray(matrix).sum(axis=0).ravel()


def batch_values(adata: ad.AnnData, batch_column: str) -> np.ndarray:
    if batch_column in adata.obs:
        return adata.obs[batch_column].astype(str).to_numpy()
    for candidate in ("batch", "measurement_set", "concat_batch"):
        if candidate in adata.obs:
            return adata.obs[candidate].astype(str).to_numpy()
    return np.repeat("all", adata.n_obs)


def barcode_set(adata: ad.AnnData, mask: np.ndarray) -> set[str]:
    return set(adata.obs_names[mask].astype(str))


def write_plot(rows: list[dict], output: Path, orientation: bool, spacer_tag: str, status: str) -> None:
    labels = [row["measurement_set"] for row in rows]
    x = np.arange(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 1.5), 5.2))
    ax.bar(x - width, [row["rna_cells"] for row in rows], width, label="RNA cells", color="#2563eb")
    ax.bar(x, [row["guide_cells"] for row in rows], width, label="guide-mapped cells", color="#f59e0b")
    ax.bar(x + width, [row["overlap_cells"] for row in rows], width, label="RNA-guide overlap", color="#10b981")
    ax.set_yscale("symlog", linthresh=1)
    ax.set_ylabel("Cells (symlog scale)")
    ax.set_xticks(x, labels, rotation=35, ha="right")
    ax.legend(frameon=False, ncol=3)
    orientation_label = "reverse-complemented" if orientation else "as supplied"
    ax.set_title(
        f"Guide mapping and barcode overlap — {status}\n"
        f"guide reference: {orientation_label}; spacer tag: {spacer_tag or 'none'}"
    )
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rna", type=Path, required=True)
    parser.add_argument("--guide", type=Path, required=True)
    parser.add_argument("--mudata", type=Path, required=True)
    parser.add_argument("--guide-metadata", type=Path, required=True)
    parser.add_argument("--reverse-complement-guides", type=parse_bool, required=True)
    parser.add_argument("--spacer-tag", default="")
    parser.add_argument("--batch-column", default="batch")
    parser.add_argument("--min-overlap-cells-per-set", type=int, default=20)
    parser.add_argument("--min-guide-to-rna-fraction", type=float, default=0.001)
    parser.add_argument("--min-overlap-to-rna-fraction", type=float, default=0.5)
    parser.add_argument("--min-overlap-to-guide-fraction", type=float, default=0.0)
    parser.add_argument("--min-recovered-guide-fraction", type=float, default=0.5)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rna = ad.read_h5ad(args.rna, backed="r")
    guide = ad.read_h5ad(args.guide)
    mdata = md.read_h5mu(args.mudata, backed="r")

    expected = read_expected_guides(args.guide_metadata)
    expected_set = set(expected)
    mapped_set = set(map(str, guide.var_names))
    feature_totals = counts_by_feature(guide.X)
    nonzero_set = set(map(str, guide.var_names[feature_totals > 0]))
    recovered_expected = expected_set & nonzero_set

    rna_batches = batch_values(rna, args.batch_column)
    guide_batches = batch_values(guide, args.batch_column)
    measurement_sets = sorted(set(rna_batches) | set(guide_batches))
    rows: list[dict] = []
    failures: list[str] = []
    for measurement_set in measurement_sets:
        rna_mask = rna_batches == measurement_set
        guide_mask = guide_batches == measurement_set
        rna_names = barcode_set(rna, rna_mask)
        guide_names = barcode_set(guide, guide_mask)
        overlap = rna_names & guide_names
        rna_cells = len(rna_names)
        guide_cells = len(guide_names)
        overlap_cells = len(overlap)
        guide_to_rna = guide_cells / rna_cells if rna_cells else 0.0
        overlap_to_guide = overlap_cells / guide_cells if guide_cells else 0.0
        overlap_to_rna = overlap_cells / rna_cells if rna_cells else 0.0
        status = "PASS"
        reasons = []
        if overlap_cells < args.min_overlap_cells_per_set:
            reasons.append(f"overlap {overlap_cells} < {args.min_overlap_cells_per_set}")
        if guide_to_rna < args.min_guide_to_rna_fraction:
            reasons.append(
                f"guide/RNA {guide_to_rna:.6f} < {args.min_guide_to_rna_fraction:.6f}"
            )
        if overlap_to_rna < args.min_overlap_to_rna_fraction:
            reasons.append(
                f"overlap/RNA {overlap_to_rna:.6f} < {args.min_overlap_to_rna_fraction:.6f}"
            )
        if overlap_to_guide < args.min_overlap_to_guide_fraction:
            reasons.append(
                f"overlap/guide {overlap_to_guide:.6f} < {args.min_overlap_to_guide_fraction:.6f}"
            )
        if reasons:
            status = "FAIL"
            failures.append(f"{measurement_set}: " + "; ".join(reasons))
        rows.append(
            {
                "measurement_set": measurement_set,
                "rna_cells": rna_cells,
                "guide_cells": guide_cells,
                "overlap_cells": overlap_cells,
                "guide_to_rna_fraction": guide_to_rna,
                "overlap_to_guide_fraction": overlap_to_guide,
                "overlap_to_rna_fraction": overlap_to_rna,
                "status": status,
                "reason": "; ".join(reasons),
            }
        )

    recovered_fraction = len(recovered_expected) / len(expected_set) if expected_set else 0.0
    if recovered_fraction < args.min_recovered_guide_fraction:
        failures.append(
            f"expected guides with nonzero counts {recovered_fraction:.6f} "
            f"< {args.min_recovered_guide_fraction:.6f}"
        )
    status = "PASS" if not failures else "FAIL"
    summary = {
        "schema_version": "1.1",
        "status": status,
        "configured_orientation": {
            "reverse_complement_guides": args.reverse_complement_guides,
            "spacer_tag": args.spacer_tag,
        },
        "thresholds": {
            "min_overlap_cells_per_set": args.min_overlap_cells_per_set,
            "min_guide_to_rna_fraction": args.min_guide_to_rna_fraction,
            "min_overlap_to_guide_fraction": args.min_overlap_to_guide_fraction,
            "min_overlap_to_rna_fraction": args.min_overlap_to_rna_fraction,
            "min_recovered_guide_fraction": args.min_recovered_guide_fraction,
        },
        "overall": {
            "rna_cells": int(rna.n_obs),
            "guide_cells": int(guide.n_obs),
            "mudata_intersection_cells": int(mdata.n_obs),
            "expected_guides": len(expected_set),
            "mapped_reference_guides": len(mapped_set),
            "expected_guides_in_reference": len(expected_set & mapped_set),
            "expected_guides_with_nonzero_counts": len(recovered_expected),
            "recovered_guide_fraction": recovered_fraction,
        },
        "measurement_sets": rows,
        "failures": failures,
    }
    (args.outdir / "guide_mapping_qc.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    with (args.outdir / "guide_mapping_qc.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    write_plot(
        rows,
        args.outdir / "guide_mapping_orientation_qc.png",
        args.reverse_complement_guides,
        args.spacer_tag,
        status,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
