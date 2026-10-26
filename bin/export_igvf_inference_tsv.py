#!/usr/bin/env python3
"""Stream the four canonical inference tables to IGVF-friendly TSV files."""

import argparse
import gzip
from pathlib import Path


NAMES = (
    "local_analysis_per_element_output",
    "local_analysis_per_guide_output",
    "global_analysis_per_element_output",
    "global_analysis_per_guide_output",
)


def convert(source: Path, destination: Path) -> None:
    if source.suffix == ".parquet":
        import polars as pl

        pl.scan_parquet(str(source)).sink_csv(str(destination), separator="\t")
    elif source.name.endswith(".tsv.gz"):
        with gzip.open(source, "rb") as reader, destination.open("wb") as writer:
            while block := reader.read(8 * 1024 * 1024):
                writer.write(block)
    elif source.suffix == ".tsv":
        with source.open("rb") as reader, destination.open("wb") as writer:
            while block := reader.read(8 * 1024 * 1024):
                writer.write(block)
    else:
        raise ValueError(f"Unsupported inference table: {source}")
    if destination.stat().st_size == 0:
        raise ValueError(f"Empty IGVF inference table: {destination}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", nargs=4, type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for expected, source in zip(NAMES, args.sources):
        if not source.name.startswith(expected + "."):
            raise ValueError(f"Expected {expected}, got {source.name}")
        convert(source, args.output_dir / f"{expected}.tsv")


if __name__ == "__main__":
    main()
