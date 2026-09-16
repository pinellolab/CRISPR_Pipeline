#!/usr/bin/env python3
"""Collect per-measurement-set RNA QC plots and audit tables."""

import argparse
import shutil
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qc_dirs", nargs="+")
    parser.add_argument("--output", type=Path, default=Path("figures"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    tables = []
    for directory in sorted(map(Path, args.qc_dirs), key=lambda path: path.name):
        for image in directory.glob("*.png"):
            shutil.copy2(image, args.output / image.name)
        for table in directory.glob("measurement_set_qc_*.tsv"):
            tables.append(pd.read_csv(table, sep="\t"))
    if not tables:
        raise ValueError("No per-measurement-set QC tables were produced")
    combined = pd.concat(tables, ignore_index=True).sort_values("measurement_set")
    combined.to_csv(args.output / "measurement_set_qc_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
