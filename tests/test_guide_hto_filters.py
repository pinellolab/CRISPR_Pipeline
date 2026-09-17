import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


ad = pytest.importorskip("anndata")
md = pytest.importorskip("mudata", exc_type=ImportError)


ROOT = Path(__file__).parents[1]


def load_script(name):
    path = ROOT / "bin" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GUIDE_QC = load_script("filter_guide_assignment_qc.py")
HTO_QC = load_script("filter_hto_post_clone.py")


def make_mudata(tmp_path):
    cells = pd.Index([f"cell{i}" for i in range(8)])
    obs = pd.DataFrame({"batch": ["A"] * 6 + ["B"] * 2}, index=cells)
    gene = ad.AnnData(sparse.csr_matrix(np.ones((8, 2))), obs=obs.copy(), var=pd.DataFrame(index=["g1", "g2"]))
    guide = ad.AnnData(sparse.csr_matrix(np.ones((8, 3))), obs=obs.copy(), var=pd.DataFrame(index=["sg1", "sg2", "sg3"]))
    guide.layers["guide_assignment"] = sparse.csr_matrix(
        [[1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 1, 0], [0, 0, 1], [1, 0, 0], [1, 0, 0], [0, 1, 0]]
    )
    hashing_obs = obs.copy()
    hashing_obs["hto_type_split"] = ["H1", "H1", "H2", "negative", "multiplets", "H2", "H3", "H3"]
    hashing = ad.AnnData(sparse.csr_matrix(np.ones((8, 3))), obs=hashing_obs, var=pd.DataFrame(index=["H1", "H2", "H3"]))
    result = md.MuData({"gene": gene, "guide": guide, "hashing": hashing})
    path = tmp_path / "input.h5mu"
    result.write_h5mu(path)
    return path


def test_guide_filter_uses_assignment_layer_and_reports_batches(tmp_path, monkeypatch):
    source = make_mudata(tmp_path)
    output = tmp_path / "guide_filtered.h5mu"
    outdir = tmp_path / "guide_qc"
    monkeypatch.setattr(
        "sys.argv",
        ["filter_guide_assignment_qc.py", str(source), str(output), "--outdir", str(outdir), "--max-guides-per-cell", "2"],
    )
    GUIDE_QC.main()
    filtered = md.read_h5mu(output)
    assert filtered.n_obs == 7
    assert "cell2" not in filtered.obs_names
    flow = pd.read_csv(outdir / "guide_assignment_filter_flow.tsv", sep="\t")
    assert flow["measurement_set"].tolist() == ["all", "A", "B"]
    assert flow.loc[flow["measurement_set"] == "A", "cells_removed"].iloc[0] == 1


def test_hto_support_is_computed_on_post_clone_input_and_keeps_singlets(tmp_path, monkeypatch):
    source = make_mudata(tmp_path)
    # Emulate upstream guide/clone filtering: only these six cells reach HTO QC.
    mdata = md.read_h5mu(source)[["cell0", "cell1", "cell2", "cell3", "cell6", "cell7"]].copy()
    post_clone = tmp_path / "post_clone.h5mu"
    mdata.write_h5mu(post_clone)
    output = tmp_path / "hto_filtered.h5mu"
    outdir = tmp_path / "hto_qc"
    monkeypatch.setattr(
        "sys.argv",
        [
            "filter_hto_post_clone.py", str(post_clone), str(output), "--outdir", str(outdir),
            "--min-positive-cells", "2", "--singlet-only", "true",
            "--filtered-hashing-output", str(tmp_path / "filtered.h5ad"),
            "--unfiltered-hashing-output", str(tmp_path / "unfiltered.h5ad"),
        ],
    )
    HTO_QC.main()
    filtered = md.read_h5mu(output)
    assert filtered.obs_names.tolist() == ["cell0", "cell1", "cell6", "cell7"]
    support = pd.read_csv(outdir / "hto_positive_cell_support.tsv", sep="\t")
    called = support.loc[support["called"], ["measurement_set", "hto_label"]]
    assert called.to_records(index=False).tolist() == [("A", "H1"), ("B", "H3")]
    metrics = (outdir / "hto_filter_metrics.json").read_text()
    assert "after_guide_assignment_and_clone_removal" in metrics
