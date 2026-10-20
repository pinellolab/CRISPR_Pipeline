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
