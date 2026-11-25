"""Fast, dependency-free contract checks for the generated documentation site."""

import csv
from collections import Counter
import json
import subprocess
import sys
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]


def test_site_sources_and_schema_reference():
    subprocess.run([sys.executable, str(ROOT / "scripts/build_docs_site.py")], cwd=ROOT, check=True)
    site = ROOT / "build/docs-site"
    assert (site / "index.md").is_file()
    assert (site / "mudata_schema.md").is_file()
    assert (site / "output-gallery.md").is_file()
    assert (site / "tool-and-module-map.md").is_file()
    for filename in (
        "garyhon_rna_knee.png",
        "huangfu_guide_mapping_overlap.png",
        "garyhon_pca_by_measurement_set.png",
        "garyhon_tf_benchmark.png",
        "huangfu_ntc_aucpr_guide.png",
    ):
        assert (site / "images/result_examples" / filename).is_file()
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

def test_current_study_examples_are_portable():
    downloads = ROOT / "build/docs-site/downloadable-examples"
    subprocess.run([sys.executable, str(ROOT / "scripts/build_docs_site.py")], cwd=ROOT, check=True)
    assert (downloads / "nextflow_cc.config").is_file()
    assert (downloads / "nextflow_tapseq.config").is_file()
    for name, expected in (
        ("gasperini_at_scale", (32, 388)),
        ("replogle_k562_essential", (48, 576)),
        ("gary_hon_tf_benchmark", (3, 12)),
    ):
        dataset = downloads / name
        with (dataset / "samplesheet.csv").open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        counts = Counter((row["measurement_sets"], row["file_modality"]) for row in rows)
        assert (len({key for key, _ in counts}), len(rows)) == expected
        for row in rows:
            for field in ("R1_path", "R2_path", "barcode_onlist", "seqspec", "guide_design", "barcode_hashtag_map"):
                value = row[field]
                assert not value or not Path(value).is_absolute()
                if value and field in ("seqspec", "guide_design", "barcode_hashtag_map"):
                    assert (dataset / value).is_file()

        with (dataset / "sample_metadata.tsv").open(newline="") as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        for sample in samples:
            for modality in ("scRNA", "gRNA", "hash"):
                assert int(sample[modality]) == counts[(sample["measurement_sets"], modality)]

        excluded = {row["barcode_onlist"] for row in rows}
        with ZipFile(downloads / f"{name}.zip") as archive:
            assert archive.testzip() is None
            members = {Path(path).name for path in archive.namelist()}
            assert {"nextflow.config", "sample_metadata.tsv", "samplesheet.csv", "prepare_samplesheet.py"} <= members
            assert not members & excluded
            assert all(".fastq" not in member and ".fq." not in member for member in members)
