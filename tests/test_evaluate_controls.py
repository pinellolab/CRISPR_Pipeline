import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from evaluate_controls import perform_binary_evaluation, run_evaluation_controls


def test_control_evaluation_skips_cleanly_without_negative_controls(tmp_path):
    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(
        X=np.ones((1, 1)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["GENE1"]),
    )
    guide_var = pd.DataFrame(
        {
            "guide_id": ["guide1"],
            "intended_target_name": ["GENE1"],
            "targeting": [True],
        },
        index=["guide1"],
    )
    guide = ad.AnnData(X=np.ones((1, 1)), obs=obs.copy(), var=guide_var)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns["global_analysis_per_guide_results"] = pd.DataFrame(
        {
            "guide_id": ["guide1"],
            "gene_id": ["GENE1"],
            "perturbo_log2_fc": [-1.0],
            "perturbo_p_value": [0.01],
        }
    )

    run_evaluation_controls(mdata, outdir=tmp_path)

    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "no non-targeting guides" in marker.read_text()


def test_binary_evaluation_skips_empty_input(tmp_path):
    assert not perform_binary_evaluation([], [], tmp_path, plot=False)
    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "valid_rows=0" in marker.read_text()


def test_control_evaluation_skips_cleanly_without_global_results(tmp_path):
    class LocalOnlyResult:
        uns = {"local_analysis_per_guide_results": pd.DataFrame()}

    run_evaluation_controls(LocalOnlyResult(), outdir=tmp_path)

    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "global PerTurbo guide results are not present" in marker.read_text()
