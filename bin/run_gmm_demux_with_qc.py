#!/usr/bin/env python3
"""Run GMM-Demux deterministically and reject degenerate HTO fits.

GMM-Demux fits one two-component Gaussian mixture per HTO.  With a weak tag,
one random initialization can separate zero from non-zero counts instead of
background from true signal.  In that failure mode every cell with a single
ambient UMI is called positive.  This wrapper tries a deterministic sequence
of seeds and accepts the first fit that passes an objective, predeclared gate.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread


def _read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def read_hto_counts(matrix_dir: Path, hto_names: list[str]) -> pd.DataFrame:
    features_path = matrix_dir / "features.tsv.gz"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    matrix_path = matrix_dir / "matrix.mtx.gz"
    features = [line.split("\t")[1] if "\t" in line else line for line in _read_lines(features_path)]
    barcodes = _read_lines(barcodes_path)
    missing = sorted(set(hto_names) - set(features))
    if missing:
        raise ValueError(f"HTOs missing from features.tsv.gz: {missing}")

    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr()
    if matrix.shape != (len(features), len(barcodes)):
        raise ValueError(
            "Matrix dimensions do not match feature/barcode files: "
            f"matrix={matrix.shape}, features={len(features)}, barcodes={len(barcodes)}"
        )
    indices = [features.index(hto) for hto in hto_names]
    return pd.DataFrame(matrix[indices, :].T.toarray(), index=barcodes, columns=hto_names)


def read_assignments(report_path: Path, config_path: Path) -> pd.DataFrame:
    report = pd.read_csv(report_path, index_col=0)
    config = pd.read_csv(config_path, header=None, names=["cluster_id", "hto_type"], skipinitialspace=True)
    config["hto_type"] = config["hto_type"].astype(str).str.strip()
    cluster_to_label = dict(zip(config["cluster_id"], config["hto_type"]))
    report["hto_type"] = report["Cluster_id"].map(cluster_to_label)
    if report["hto_type"].isna().any():
        missing = sorted(report.loc[report["hto_type"].isna(), "Cluster_id"].unique().tolist())
        raise ValueError(f"GMM report references cluster IDs absent from config: {missing}")
    return report


def evaluate_fit(
    counts: pd.DataFrame,
    assignments: pd.DataFrame,
    reject_nonzero_positive: bool = True,
    equivalence_fraction: float = 0.995,
) -> dict:
    missing = assignments.index.difference(counts.index)
    if len(missing):
        raise ValueError(f"GMM report contains {len(missing)} barcodes absent from the input matrix")
    aligned = counts.loc[assignments.index]
    labels = assignments["hto_type"].astype(str).str.split("-")
    per_hto = {}
    rejected = []

    for hto in counts.columns:
        values = aligned[hto].to_numpy()
        positive = labels.apply(lambda parts: hto in parts).to_numpy(dtype=bool)
        nonzero = values > 0
        nonzero_n = int(nonzero.sum())
        zero_n = int((~nonzero).sum())
        positive_nonzero_fraction = float((positive & nonzero).sum() / nonzero_n) if nonzero_n else 0.0
        positive_zero_fraction = float((positive & ~nonzero).sum() / zero_n) if zero_n else 0.0
        positive_values = values[positive]
        min_positive_count = float(positive_values.min()) if positive_values.size else None
        exact_nonzero_equivalence = bool(np.array_equal(positive, nonzero))
        near_nonzero_equivalence = bool(
            positive_values.size
            and min_positive_count <= 1
            and positive_nonzero_fraction >= equivalence_fraction
            and positive_zero_fraction <= (1.0 - equivalence_fraction)
        )
        reason = None
        if reject_nonzero_positive and near_nonzero_equivalence:
            reason = (
                f"{hto}: positive calls reproduce count>0 "
                f"(nonzero-positive={positive_nonzero_fraction:.4f}, "
                f"zero-positive={positive_zero_fraction:.4f}, "
                f"minimum-positive-count={min_positive_count:g})"
            )
            rejected.append(reason)
        per_hto[hto] = {
            "positive_cells": int(positive.sum()),
            "nonzero_cells": nonzero_n,
            "minimum_positive_count": min_positive_count,
            "positive_nonzero_fraction": positive_nonzero_fraction,
            "positive_zero_fraction": positive_zero_fraction,
            "exact_nonzero_equivalence": exact_nonzero_equivalence,
            "near_nonzero_equivalence": near_nonzero_equivalence,
            "rejection_reason": reason,
        }

    return {"accepted": not rejected, "rejection_reasons": rejected, "per_hto": per_hto}


def run_attempt(
    executable: str,
    matrix_dir: Path,
    hto_string: str,
    seed: int,
    attempt_dir: Path,
    ssd_dir: Path,
) -> None:
    shutil.rmtree(attempt_dir, ignore_errors=True)
    shutil.rmtree(ssd_dir, ignore_errors=True)
    subprocess.run(
        [
            executable,
            str(matrix_dir),
            hto_string,
            "-f",
            str(attempt_dir),
            "-o",
            str(ssd_dir),
            "--random_seed",
            str(seed),
        ],
        check=True,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir", required=True, type=Path)
    parser.add_argument("--hto-names", required=True, help="Comma-separated HTO names in matrix order")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--ssd-output-dir", required=True, type=Path)
    parser.add_argument("--qc-json", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-attempts", type=int, default=10)
    parser.add_argument("--gmm-executable", default="GMM-demux")
    parser.add_argument("--reject-nonzero-positive", action="store_true")
    args = parser.parse_args()

    if args.max_attempts < 1:
        parser.error("--max-attempts must be at least 1")
    hto_names = [value.strip() for value in args.hto_names.split(",") if value.strip()]
    if not hto_names:
        parser.error("--hto-names must contain at least one HTO")

    counts = read_hto_counts(args.matrix_dir, hto_names)
    attempts = []
    selected = None
    for offset in range(args.max_attempts):
        seed = args.seed + offset
        attempt_dir = Path(f"GMM_ATTEMPT_{offset:02d}_SEED_{seed}")
        ssd_dir = Path(f"SSD_ATTEMPT_{offset:02d}_SEED_{seed}")
        run_attempt(args.gmm_executable, args.matrix_dir, args.hto_names, seed, attempt_dir, ssd_dir)
        assignments = read_assignments(attempt_dir / "GMM_full.csv", attempt_dir / "GMM_full.config")
        assessment = evaluate_fit(counts, assignments, args.reject_nonzero_positive)
        assessment.update({"attempt": offset + 1, "seed": seed})
        attempts.append(assessment)
        if assessment["accepted"]:
            selected = (seed, attempt_dir, ssd_dir)
            break

    payload = {
        "base_seed": args.seed,
        "max_attempts": args.max_attempts,
        "reject_nonzero_positive": args.reject_nonzero_positive,
        "selected_seed": selected[0] if selected else None,
        "attempts": attempts,
    }
    args.qc_json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if selected is None:
        reasons = "; ".join(reason for attempt in attempts for reason in attempt["rejection_reasons"])
        raise RuntimeError(
            f"No acceptable GMM-Demux fit after {args.max_attempts} deterministic attempts. {reasons}"
        )

    _, selected_dir, selected_ssd = selected
    shutil.rmtree(args.output_dir, ignore_errors=True)
    shutil.rmtree(args.ssd_output_dir, ignore_errors=True)
    shutil.copytree(selected_dir, args.output_dir)
    shutil.copytree(selected_ssd, args.ssd_output_dir)
    print(f"Selected deterministic GMM-Demux seed {selected[0]} after {len(attempts)} attempt(s)")


if __name__ == "__main__":
    main()
