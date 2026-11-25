#!/usr/bin/env python3
"""Resolve the paths in a versioned CRISPR Pipeline example samplesheet."""

import argparse
import csv
from collections import Counter
from pathlib import Path


def prepare(template, fastq_dir, barcode_dir, output):
    template = template.resolve(strict=True)
    output = output.resolve()
    if output == template:
        raise ValueError("Choose a new path for the resolved samplesheet.")

    with template.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames
        rows = list(reader)
    if not rows:
        raise ValueError(f"No rows in {template}")

    roots = {
        "R1_path": fastq_dir, "R2_path": fastq_dir, "barcode_onlist": barcode_dir,
        "seqspec": template.parent, "guide_design": template.parent,
        "barcode_hashtag_map": template.parent,
    }
    missing = set()
    for row in rows:
        for field, root in roots.items():
            if not row.get(field):
                continue
            path = (root / row[field]).absolute()
            row[field] = str(path)
            if not path.is_file():
                missing.add(path)
    if missing:
        shown = "\n".join(f"  {path}" for path in sorted(missing)[:10])
        raise ValueError(f"{len(missing)} input files are missing. First missing paths:\n{shown}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = Counter(row["file_modality"] for row in rows)
    print(f"Wrote {output} with {len(rows)} rows ({dict(counts)})")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("samplesheet", type=Path, help="Example samplesheet.csv")
    parser.add_argument("--fastq-dir", type=Path, required=True)
    parser.add_argument("--barcode-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        prepare(args.samplesheet, args.fastq_dir, args.barcode_dir, args.output)
    except (OSError, ValueError) as error:
        parser.exit(2, f"{error}\n")


if __name__ == "__main__":
    main()
