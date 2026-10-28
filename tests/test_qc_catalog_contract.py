"""Keep exported definitions and emitted guide QC statistics synchronized."""
import csv
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "bin"))
from export_qc_metric_catalog import flatten_catalog, write_json, write_markdown, write_tsv
from qc_metrics_json import METRIC_CATALOG


def test_generated_catalogs_match_committed_files(tmp_path):
    records = flatten_catalog()
    for filename, writer in [("qc_outputs_flat.tsv", write_tsv),
                             ("qc_outputs_flat.json", write_json),
                             ("qc_outputs.md", write_markdown)]:
        generated = tmp_path / filename
        writer(generated, records)
        assert generated.read_bytes() == (ROOT / "docs" / filename).read_bytes(), filename


def test_portal_mapping_and_execution_status_fields_are_documented():
    keys = {r["flat_key"] for r in flatten_catalog()}
    expected = {
        "guide_mapping_qc.overlap_to_rna_fraction",
        "additional_qc_guide.n_cells_before_assigned_guide_filter",
        "additional_qc_guide.n_cells_without_assigned_guide",
        "additional_qc_guide.assigned_guide_filter_applied",
        "additional_qc_guide.assigned_guide_counts_source",
        "measurement_set_rna_qc.scrublet_status",
        "measurement_set_rna_qc.scrublet_skip_reason",
        "measurement_set_rna_qc.scrublet_failure_policy",
    }
    assert expected <= keys


def test_guide_emitter_supplies_catalog_statistics():
    ad = pytest.importorskip("anndata")
    pytest.importorskip("mudata")
    from mapping_guide import compute_guide_metrics

    calls = np.zeros((8, 4))
    for col, count in enumerate([0, 1, 3, 8]):
        calls[:count, col] = 1
    guide = ad.AnnData(calls)
    guide.obs["total_guide_umis"] = calls.sum(axis=1)
    guide.obs["n_guides_per_cell"] = calls.sum(axis=1)
    guide.var["n_cells_per_guide"] = calls.sum(axis=0)
    row = compute_guide_metrics(guide)
    declared = {m["name"] for m in METRIC_CATALOG["additional_qc_guide"]["metrics"]}
    guide_stats = {k for k in declared if k.startswith("cells_per_guide_")}
    assert guide_stats <= row.keys()
    assert row["cells_per_guide_q25"] == pytest.approx(0.75)
    assert row["cells_per_guide_q75"] == pytest.approx(4.25)
    batch_row = compute_guide_metrics(guide, batch_label="set1", include_per_guide_stats=False)
    assert "cells_per_guide_q25" not in batch_row
    assert "cells_per_guide_q75" not in batch_row


def test_mapping_emitter_columns_match_catalog_and_overlap_denominator(tmp_path, monkeypatch):
    ad = pytest.importorskip("anndata")
    md = pytest.importorskip("mudata")
    import guide_mapping_qc

    rna = ad.AnnData(np.ones((3, 1)))
    rna.obs_names = ["a", "b", "c"]
    guide = ad.AnnData(np.ones((2, 1)))
    guide.obs_names = ["b", "d"]
    guide.var_names = ["g1"]
    rna.write_h5ad(tmp_path / "rna.h5ad")
    guide.write_h5ad(tmp_path / "guide.h5ad")
    md.MuData({"gene": rna[["b"], :].copy(), "guide": guide[["b"], :].copy()}).write_h5mu(tmp_path / "input.h5mu")
    (tmp_path / "guides.tsv").write_text("guide_id\ng1\n")
    monkeypatch.setattr(guide_mapping_qc, "write_plot", lambda *args: None)
    monkeypatch.setattr(sys, "argv", ["guide_mapping_qc", "--rna", str(tmp_path / "rna.h5ad"),
        "--guide", str(tmp_path / "guide.h5ad"), "--mudata", str(tmp_path / "input.h5mu"),
        "--guide-metadata", str(tmp_path / "guides.tsv"), "--reverse-complement-guides", "false",
        "--min-overlap-cells-per-set", "0", "--min-overlap-to-rna-fraction", "0",
        "--outdir", str(tmp_path / "qc")])
    assert guide_mapping_qc.main() == 0
    with (tmp_path / "qc/guide_mapping_qc.tsv").open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    declared = {m["name"] for m in METRIC_CATALOG["guide_mapping_qc"]["metrics"]}
    assert set(row) == declared
    assert float(row["overlap_to_rna_fraction"]) == pytest.approx(1 / 3)
    assert float(row["overlap_to_guide_fraction"]) == pytest.approx(1 / 2)
