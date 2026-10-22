"""Fast, dependency-free contract checks for the generated documentation site."""

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_site_sources_and_schema_reference():
    subprocess.run([sys.executable, str(ROOT / "scripts/build_docs_site.py")], cwd=ROOT, check=True)
    site = ROOT / "build/docs-site"
    assert (site / "index.md").is_file()
    assert (site / "mudata_schema.md").is_file()
    assert (site / "downloadable-examples/hashing_samplesheet.csv").is_file()
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    expected = {name for section in schema["$defs"].values() for name in section.get("properties", {})}
    reference = (site / "parameters.md").read_text()
    assert len(expected) >= 100
    assert all(f"| `{name}` |" in reference for name in expected)


def test_examples_follow_current_samplesheet_header():
    required = {"R1_path", "file_modality", "measurement_sets", "sequencing_run", "lane", "seqspec", "barcode_onlist", "guide_design"}
    for filename in ("minimal_samplesheet.csv", "hashing_samplesheet.csv"):
        with (ROOT / "docs/examples" / filename).open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert rows
        assert required.issubset(rows[0])
        assert all(row["R1_path"] and row["barcode_onlist"] and row["guide_design"] for row in rows)
        modalities = {row["file_modality"] for row in rows}
        assert {"scRNA", "gRNA"}.issubset(modalities)
        if "hash" in modalities:
            assert all(row["barcode_hashtag_map"] for row in rows)
