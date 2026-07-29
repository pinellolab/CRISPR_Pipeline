#!/usr/bin/env python
"""Export the pipeline QC metric catalog as flat JSON, TSV, and Markdown."""

import argparse
import csv
import json
from pathlib import Path

from qc_metrics_json import METRIC_CATALOG


def flatten_catalog():
    records = []
    for section, section_spec in METRIC_CATALOG.items():
        source = section_spec.get("source_artifact")
        if source is None:
            source = section_spec.get("source_artifact_pattern")
        row_level = section_spec.get("row_level", "measurement-set or run-level")
        for metric in section_spec["metrics"]:
            records.append(
                {
                    "flat_key": f"{section}.{metric['name']}",
                    "section": section,
                    "metric": metric["name"],
                    "description": metric["description"],
                    "unit": metric.get("unit"),
                    "row_level": row_level,
                    "source_artifact": source,
                }
            )
    return records


def write_json(path, records):
    path.write_text(json.dumps(records, indent=2) + "\n")


def write_tsv(path, records):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(records[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)


def markdown_value(value):
    if value is None:
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ")


def write_markdown(path, records):
    lines = [
        "# Pipeline QC outputs",
        "",
        (
            "This table is generated from `bin/qc_metrics_json.py`. The flat key "
            "combines the catalog section and metric name; it remains unique even "
            "when different QC artifacts use the same column name."
        ),
        "",
        "| Flat key | Description | Unit | Row level | Source artifact |",
        "|---|---|---|---|---|",
    ]
    for record in records:
        lines.append(
            "| {flat_key} | {description} | {unit} | {row_level} | "
            "{source_artifact} |".format(
                **{key: markdown_value(value) for key, value in record.items()}
            )
        )
    lines.append("")
    path.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records = flatten_catalog()
    write_json(args.output_dir / "qc_outputs_flat.json", records)
    write_tsv(args.output_dir / "qc_outputs_flat.tsv", records)
    write_markdown(args.output_dir / "qc_outputs.md", records)
    print(f"Wrote {len(records)} QC metric definitions to {args.output_dir}")


if __name__ == "__main__":
    main()
