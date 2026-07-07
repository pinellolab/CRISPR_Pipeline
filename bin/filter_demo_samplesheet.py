#!/usr/bin/env python3

"""Create a pre-run samplesheet containing one complete measurement set."""

import argparse
import csv
import sys
from pathlib import Path


SUPPORTED_MODALITIES = {"scrna", "grna", "hash"}


def detect_delimiter(path: Path) -> str:
    with path.open(newline="") as handle:
        header = next((line for line in handle if line.strip()), "")
    return "\t" if "\t" in header else ","


def count_scrna_fastq_files(rows) -> int:
    return len(
        {
            str(row.get(column)).strip()
            for row in rows
            if str(row.get("file_modality") or "").strip().lower() == "scrna"
            for column in ("R1_path", "R2_path")
            if str(row.get(column) or "").strip()
        }
    )


def filter_demo_samplesheet(input_path, output_path, require_hash=False):
    input_path = Path(input_path)
    output_path = Path(output_path)
    delimiter = detect_delimiter(input_path)

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError("Input samplesheet has no header.")
        missing = {"file_modality", "measurement_sets"}.difference(reader.fieldnames)
        if missing:
            raise ValueError(
                "Input samplesheet is missing required demo-mode column(s): "
                + ", ".join(sorted(missing))
            )
        rows = [row for row in reader if any(str(value or "").strip() for value in row.values())]
        fieldnames = reader.fieldnames

    required_modalities = {"scrna", "grna"}
    if require_hash:
        required_modalities.add("hash")

    measurement_set_order = []
    modalities_by_set = {}
    for row in rows:
        measurement_set = str(row.get("measurement_sets") or "").strip()
        modality = str(row.get("file_modality") or "").strip().lower()
        if not measurement_set:
            continue
        if measurement_set not in modalities_by_set:
            measurement_set_order.append(measurement_set)
            modalities_by_set[measurement_set] = set()
        if modality in SUPPORTED_MODALITIES:
            modalities_by_set[measurement_set].add(modality)

    complete_sets = [
        measurement_set
        for measurement_set in measurement_set_order
        if required_modalities.issubset(modalities_by_set[measurement_set])
    ]
    rows_by_set = {
        measurement_set: [
            row
            for row in rows
            if str(row.get("measurement_sets") or "").strip() == measurement_set
        ]
        for measurement_set in complete_sets
    }
    selected_set = (
        min(
            complete_sets,
            key=lambda measurement_set: (
                count_scrna_fastq_files(rows_by_set[measurement_set]),
                measurement_set,
            ),
        )
        if complete_sets
        else None
    )
    if selected_set is None:
        required = ", ".join(sorted(required_modalities))
        observed = "; ".join(
            f"{measurement_set}: {', '.join(sorted(modalities_by_set[measurement_set])) or '<none>'}"
            for measurement_set in measurement_set_order
        )
        raise ValueError(
            "Demo mode could not find a complete measurement set with modalities "
            f"[{required}]. Observed [{observed or '<no measurement sets>'}]."
        )

    selected_rows = [
        row
        for row in rows
        if str(row.get("measurement_sets") or "").strip() == selected_set
        and str(row.get("file_modality") or "").strip().lower() in SUPPORTED_MODALITIES
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter=delimiter,
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(selected_rows)

    return selected_set, selected_rows


def write_warning(path, selected_set, row_count, scrna_fastq_file_count):
    Path(path).write_text(
        "DEMO MODE / PRE-RUN ONLY\n"
        "========================\n"
        "This pipeline run used a filtered samplesheet and is not a final-results run.\n"
        f"Selected measurement set: {selected_set}\n"
        f"scRNA FASTQ files retained: {scrna_fastq_file_count}\n"
        f"Samplesheet rows retained: {row_count}\n"
        "Do not report, publish, or treat these outputs as results from the full dataset.\n"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Keep the smallest complete RNA/guide measurement set for a demo pre-run."
    )
    parser.add_argument("--input", required=True, help="Input CSV or TSV samplesheet")
    parser.add_argument("--output", required=True, help="Filtered CSV or TSV samplesheet")
    parser.add_argument(
        "--require-hash",
        action="store_true",
        help="Require the selected measurement set to contain hash rows",
    )
    parser.add_argument(
        "--warning-output",
        default="DEMO_MODE_WARNING.txt",
        help="Path for the prominent pre-run warning marker",
    )
    args = parser.parse_args()

    try:
        selected_set, selected_rows = filter_demo_samplesheet(
            args.input,
            args.output,
            require_hash=args.require_hash,
        )
    except ValueError as exc:
        parser.error(str(exc))

    scrna_fastq_file_count = count_scrna_fastq_files(selected_rows)
    write_warning(
        args.warning_output,
        selected_set,
        len(selected_rows),
        scrna_fastq_file_count,
    )
    print(
        "WARNING: DEMO MODE retained only measurement set "
        f"'{selected_set}' ({scrna_fastq_file_count} scRNA FASTQ files; "
        f"{len(selected_rows)} rows). PRE-RUN ONLY; not final results.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
