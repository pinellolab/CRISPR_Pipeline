import csv
import sys
from pathlib import Path

import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from filter_demo_samplesheet import filter_demo_samplesheet, write_warning


FIELDNAMES = ["R1_path", "R2_path", "file_modality", "measurement_sets"]


def write_samplesheet(path, rows, delimiter=","):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def test_selects_complete_set_with_fewest_scrna_fastq_files(tmp_path):
    source = tmp_path / "samples.csv"
    output = tmp_path / "demo.csv"
    rows = [
        {"R1_path": "b-rna", "file_modality": "scRNA", "measurement_sets": "B"},
        {"R1_path": "b-guide", "file_modality": "gRNA", "measurement_sets": "B"},
        {"R1_path": "a-rna-r1", "R2_path": "a-rna-r2", "file_modality": "scRNA", "measurement_sets": "A"},
        {"R1_path": "a-guide", "file_modality": "gRNA", "measurement_sets": "A"},
        {"R1_path": "a-other", "file_modality": "ATAC", "measurement_sets": "A"},
    ]
    write_samplesheet(source, rows)

    selected, selected_rows = filter_demo_samplesheet(source, output)

    assert selected == "B"
    assert [row["R1_path"] for row in selected_rows] == ["b-rna", "b-guide"]
    with output.open(newline="") as handle:
        assert len(list(csv.DictReader(handle))) == 2


def test_breaks_scrna_file_count_ties_by_measurement_set_id(tmp_path):
    source = tmp_path / "samples.csv"
    output = tmp_path / "demo.csv"
    rows = [
        {"R1_path": "b-rna", "file_modality": "scRNA", "measurement_sets": "B"},
        {"R1_path": "b-guide", "file_modality": "gRNA", "measurement_sets": "B"},
        {"R1_path": "a-rna", "file_modality": "scRNA", "measurement_sets": "A"},
        {"R1_path": "a-guide", "file_modality": "gRNA", "measurement_sets": "A"},
    ]
    write_samplesheet(source, rows)

    selected, _ = filter_demo_samplesheet(source, output)

    assert selected == "A"


def test_skips_incomplete_set_and_requires_hash_when_enabled(tmp_path):
    source = tmp_path / "samples.tsv"
    output = tmp_path / "demo.tsv"
    rows = [
        {"R1_path": "a-rna", "file_modality": "scRNA", "measurement_sets": "A"},
        {"R1_path": "a-guide", "file_modality": "gRNA", "measurement_sets": "A"},
        {"R1_path": "b-rna", "file_modality": "scRNA", "measurement_sets": "B"},
        {"R1_path": "b-guide", "file_modality": "gRNA", "measurement_sets": "B"},
        {"R1_path": "b-hash", "file_modality": "hash", "measurement_sets": "B"},
    ]
    write_samplesheet(source, rows, delimiter="\t")

    selected, selected_rows = filter_demo_samplesheet(source, output, require_hash=True)

    assert selected == "B"
    assert {row["file_modality"] for row in selected_rows} == {"scRNA", "gRNA", "hash"}
    assert "\t" in output.read_text().splitlines()[0]


def test_fails_when_no_complete_measurement_set_exists(tmp_path):
    source = tmp_path / "samples.csv"
    rows = [
        {"R1_path": "a-rna", "file_modality": "scRNA", "measurement_sets": "A"},
        {"R1_path": "b-guide", "file_modality": "gRNA", "measurement_sets": "B"},
    ]
    write_samplesheet(source, rows)

    with pytest.raises(ValueError, match="could not find a complete measurement set"):
        filter_demo_samplesheet(source, tmp_path / "demo.csv")


def test_warning_marks_outputs_as_not_final(tmp_path):
    warning = tmp_path / "DEMO_MODE_WARNING.txt"

    write_warning(warning, "A", 3, 2)

    text = warning.read_text()
    assert "PRE-RUN ONLY" in text
    assert "not a final-results run" in text
    assert "Selected measurement set: A" in text
    assert "scRNA FASTQ files retained: 2" in text
