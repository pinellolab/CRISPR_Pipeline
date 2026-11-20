import json
import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

BIN = pathlib.Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))

import preinference_checkpoint
import prepare_external_mudata
from inference_covariates import declared_external_covariates, materialize_shared_covariates


def make_mudata(path):
    cells = [f"cell{i}" for i in range(4)]
    obs = pd.DataFrame(index=cells)
    gene = ad.AnnData(
        X=sparse.csr_matrix([[3, 0], [2, 1], [1, 2], [0, 3]]),
        obs=obs.copy(), var=pd.DataFrame(index=["gene1", "gene2"]),
    )
    guide = ad.AnnData(
        X=sparse.csr_matrix([[2, 0], [1, 0], [0, 2], [0, 1]]),
        obs=obs.copy(), var=pd.DataFrame(index=["guide1", "guide2"]),
    )
    guide.layers["guide_assignment"] = sparse.csr_matrix([[1, 0], [1, 0], [0, 1], [0, 1]])
    mu.MuData({"gene": gene, "guide": guide}).write_h5mu(path)


def invoke(monkeypatch, function, *args):
    monkeypatch.setattr(sys, "argv", [function.__module__, *map(str, args)])
    function()


def test_checkpoint_external_filter_and_shared_covariate(tmp_path, monkeypatch):
    source = tmp_path / "source.h5mu"
    make_mudata(source)
    gtf = tmp_path / "reference.gtf"
    gtf.write_text("test reference\n")
    checkpoint = tmp_path / "preinference_mudata.h5mu"
    manifest = tmp_path / "preinference_manifest.json"
    invoke(monkeypatch, preinference_checkpoint.main,
           "--input", source, "--gtf", gtf, "--output", checkpoint,
           "--manifest", manifest, "--run-name", "pilot", "--commit", "test")
    assert json.loads(manifest.read_text())["status"] == "READY_FOR_INFERENCE"

    keep = tmp_path / "keep.tsv"
    keep.write_text("cell_id\ncell0\ncell1\ncell2\n")
    covariates = tmp_path / "covariates.tsv"
    covariates.write_text("cell_id\texternal_quality\ncell0\t1\ncell1\t2\ncell2\t3\n")
    output = tmp_path / "analysis_input.h5mu"
    report = tmp_path / "external_mudata_report.json"
    invoke(monkeypatch, prepare_external_mudata.main,
           "--parent-manifest", manifest, "--parent-mudata", checkpoint,
           "--gtf", gtf, "--keep-cells", keep, "--covariates", covariates,
           "--output", output, "--report", report)
    result = mu.read_h5mu(output)
    assert result["gene"].obs_names.tolist() == ["cell0", "cell1", "cell2"]
    assert result["guide"].obs_names.tolist() == result["gene"].obs_names.tolist()
    assert declared_external_covariates(result) == ["external_quality"]
    assert "external_quality" in materialize_shared_covariates(result)
    assert result.obs["external_quality"].tolist() == [1.0, 2.0, 3.0]
    assert json.loads(report.read_text())["accepted_cells"] == 3


def test_reject_changed_raw_counts(tmp_path, monkeypatch):
    source = tmp_path / "source.h5mu"
    make_mudata(source)
    gtf = tmp_path / "reference.gtf"
    gtf.write_text("test reference\n")
    checkpoint = tmp_path / "preinference_mudata.h5mu"
    manifest = tmp_path / "preinference_manifest.json"
    invoke(monkeypatch, preinference_checkpoint.main,
           "--input", source, "--gtf", gtf, "--output", checkpoint,
           "--manifest", manifest, "--run-name", "pilot", "--commit", "test")
    edited = mu.read_h5mu(checkpoint)
    edited["gene"].X[0, 0] = 99
    changed = tmp_path / "changed.h5mu"
    edited.write_h5mu(changed)
    with pytest.raises(ValueError, match="raw counts"):
        invoke(monkeypatch, prepare_external_mudata.main,
               "--parent-manifest", manifest, "--parent-mudata", checkpoint,
               "--input-mudata", changed, "--gtf", gtf,
               "--output", tmp_path / "analysis_input.h5mu")


def test_embedded_external_covariate_without_tsv(tmp_path, monkeypatch):
    source = tmp_path / "source.h5mu"
    make_mudata(source)
    gtf = tmp_path / "reference.gtf"
    gtf.write_text("test reference\n")
    checkpoint = tmp_path / "preinference_mudata.h5mu"
    manifest = tmp_path / "preinference_manifest.json"
    invoke(monkeypatch, preinference_checkpoint.main,
           "--input", source, "--gtf", gtf, "--output", checkpoint,
           "--manifest", manifest, "--run-name", "pilot", "--commit", "test")

    edited = mu.read_h5mu(checkpoint)
    edited["gene"].obs["external_quality"] = [1.0, 2.0, 3.0, 4.0]
    edited.uns["external_covariates"] = ["external_quality"]
    edited_path = tmp_path / "edited.h5mu"
    edited.write_h5mu(edited_path)
    output = tmp_path / "analysis_input.h5mu"
    report = tmp_path / "external_mudata_report.json"
    invoke(monkeypatch, prepare_external_mudata.main,
           "--parent-manifest", manifest, "--parent-mudata", checkpoint,
           "--input-mudata", edited_path, "--gtf", gtf,
           "--output", output, "--report", report)
    result = mu.read_h5mu(output)
    assert declared_external_covariates(result) == ["external_quality"]
    assert result.obs["external_quality"].tolist() == [1.0, 2.0, 3.0, 4.0]
    assert "external_quality" in materialize_shared_covariates(result)
    assert json.loads(report.read_text())["external_covariates"] == ["external_quality"]


def test_continued_dashboard_keeps_parent_qc_and_shows_resumed_results(tmp_path, monkeypatch):
    import checkpoint_dashboard
    import wandb_html_monitor

    parent = tmp_path / "parent" / "preinference_checkpoint" / "checkpoint-a"
    parent.mkdir(parents=True)
    source = tmp_path / "source.h5mu"
    make_mudata(source)
    gtf = tmp_path / "reference.gtf"
    gtf.write_text("reference\n")
    manifest = parent / "preinference_manifest.json"
    checkpoint = parent / "preinference_mudata.h5mu"
    invoke(monkeypatch, preinference_checkpoint.main, "--input", source, "--gtf", gtf,
           "--output", checkpoint, "--manifest", manifest, "--run-name", "checkpoint-a", "--commit", "test")
    original_manifest = manifest.read_bytes()
    metrics = tmp_path / "pipeline_qc_metrics.json"
    metrics.write_text(json.dumps({
        "observed_metrics": {
            "mapping_json": [{"measurement_set": "Parent library", "modality": "trans",
                              "metrics": {"n_processed": 43821}}],
            "measurement_set_rna_qc": {"rows": [{"measurement_set": "Parent library", "retained_cells": 4}]},
        },
        "dashboard_summary": {"filtering_summary": {"final_mudata": {"cells": 4}}},
    }))
    seqspec = tmp_path / "guide_position_table.csv"
    seqspec.write_text("Sample,Config,TotalHits\nParent library,config-a,42\n")
    figures = tmp_path / "figures"
    figures.mkdir()
    (figures / "rna_qc_full_panel.png").write_bytes(b"parent QC plot")
    (figures / "unrelated.h5mu").write_bytes(b"must not be copied")
    invoke(monkeypatch, checkpoint_dashboard.main, "snapshot", "--manifest", manifest,
           "--metrics", metrics, "--seqspec", seqspec, "--output", parent / "checkpoint_qc", figures)
    assert not (parent / "checkpoint_qc" / "figures" / "unrelated.h5mu").exists()

    continued = tmp_path / "continued"
    continued.mkdir()
    keep = tmp_path / "keep.tsv"
    keep.write_text("cell_id\ncell0\ncell1\ncell2\n")
    accepted = continued / "analysis_input.h5mu"
    report = continued / "external_mudata_report.json"
    invoke(monkeypatch, prepare_external_mudata.main, "--parent-manifest", manifest,
           "--parent-mudata", checkpoint, "--gtf", gtf, "--keep-cells", keep,
           "--output", accepted, "--report", report)
    invoke(monkeypatch, checkpoint_dashboard.main, "restore", "--report", report, "--source", parent,
           "--guide-report", parent / "missing.json", "--output", continued / "checkpoint_dashboard")
    assert manifest.read_bytes() == original_manifest
    assert (continued / "checkpoint_dashboard" / metrics.name).read_bytes() == metrics.read_bytes()

    trace = continued / "trace.tsv"
    trace.write_text("process\tstatus\nprepare_external_mudata\tCOMPLETED\ninference_pipeline:inference_sceptre\tRUNNING\n")
    args = type("Args", (), {
        "outdir": continued, "run_name": "continuation-a", "trace": trace,
        "source_run_id": "continued", "dashboard_html": continued / "live.html",
        "nextflow_log": None, "tail_lines": 20, "max_image_bytes": 2_000_000,
    })()
    wandb_html_monitor.render_snapshot(args, "running", False)
    live = args.dashboard_html.read_text()
    assert "Continuation from checkpoint-a" in live
    assert "Checkpoint cells" in live and "Inference input cells" in live
    assert "43,821" in live and "rna_qc_full_panel" in live
    assert 'class="node checkpoint" data-family="mapping"' in live
    assert 'class="node running" data-family="inference"' in live

    results = mu.read_h5mu(accepted)
    results.uns["local_analysis_per_guide_results"] = {
        "gene_name": np.array(["first_gene", "lowest_pvalue_gene"]),
        "sceptre_p_value": np.array([0.2, 0.01]),
    }
    results.uns["global_analysis_per_guide_results"] = {
        "gene_name": np.array(["global_gene"]),
        "perturbo_p_value": np.array([0.03]),
    }
    results_path = continued / "inference_mudata.h5mu"
    results.write_h5mu(results_path)
    evaluation = continued / "evaluation_output"
    evaluation.mkdir()
    (evaluation / "evaluation_skipped.txt").write_text("Synthetic evaluation output")
    invoke(monkeypatch, checkpoint_dashboard.main, "complete", "--root", continued,
           "--mudata", results_path, "--run-name", "continuation-a")
    published = continued / "pipeline_dashboard"
    (continued / "dashboard.html").rename(published / "dashboard.html")
    (continued / "inference_preview.html").rename(published / "inference_preview.html")
    standalone = (published / "dashboard.html").read_text()
    assert "Synthetic evaluation output" in standalone
    assert "43,821" in standalone and "rna_qc_full_panel" in standalone
    assert standalone.index("lowest_pvalue_gene") < standalone.index("first_gene")
    assert "global_gene" in standalone
    trace.write_text(trace.read_text().replace("RUNNING", "COMPLETED"))
    wandb_html_monitor.render_snapshot(args, "completed", True)
    final = (published / "dashboard.html").read_text()
    assert "43,821" in final and "lowest_pvalue_gene" in final
    assert not (published / "legacy_dashboard.html").exists()

    args.outdir = parent.parents[1]
    args.run_name = "checkpoint-a"
    args.dashboard_html = args.outdir / "live.html"
    wandb_html_monitor.render_snapshot(args, "completed", True)
    assert "READY FOR INFERENCE" in args.dashboard_html.read_text()
    assert "awaiting continuation" in args.dashboard_html.read_text()


def test_older_checkpoint_retains_available_qc_and_explicit_missing_state(tmp_path, monkeypatch):
    import checkpoint_dashboard
    import wandb_html_monitor

    report = tmp_path / "external_mudata_report.json"
    report.write_text(json.dumps({
        "parent_run_name": "old", "parent_sha256": "abc", "parent_cells": 4,
        "submitted_cells": 4, "accepted_cells": 3,
    }))
    for with_qc in (False, True):
        source = tmp_path / f"source-{with_qc}"
        source.mkdir()
        if with_qc:
            qc = source / "measurement_set_qc"
            qc.mkdir()
            (qc / "measurement_set_qc_metrics.tsv").write_text(
                "measurement_set\tretained_cells\nLegacy library\t4\n"
            )
        child = tmp_path / f"child-{with_qc}"
        child.mkdir()
        invoke(monkeypatch, checkpoint_dashboard.main, "restore", "--report", report,
               "--source", source, "--guide-report", source / "missing.json",
               "--output", child / "checkpoint_dashboard")
        args = type("Args", (), {
            "outdir": child, "run_name": "continuation", "trace": child / "trace.tsv",
            "source_run_id": "continued", "dashboard_html": child / "live.html",
            "nextflow_log": None, "tail_lines": 20, "max_image_bytes": 2_000_000,
        })()
        wandb_html_monitor.render_snapshot(args, "running", False)
        content = args.dashboard_html.read_text()
        assert "Inference input cells" in content
        if with_qc:
            assert "partial published QC" in content and "Legacy library" in content
        else:
            assert "Published mapping and QC files for this checkpoint are unavailable." in content
