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


def test_control_evaluation_matches_symbol_target_to_ensembl_result(tmp_path, monkeypatch):
    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(X=np.ones((1, 1)), obs=obs.copy(), var=pd.DataFrame(index=["ENSG1"]))
    guide_var = pd.DataFrame(
        {
            "guide_id": ["targeting-guide", "control-guide"],
            "intended_target_name": ["GENE1", "non-targeting"],
            "targeting": [True, False],
        },
        index=["targeting-guide", "control-guide"],
    )
    guide = ad.AnnData(X=np.ones((1, 2)), obs=obs.copy(), var=guide_var)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns["global_analysis_per_guide_results"] = pd.DataFrame(
        {
            "guide_id": ["targeting-guide", "control-guide"],
            "gene_id": ["ENSG1", "ENSG1"],
            "gene_name": ["GENE1", "GENE1"],
            "perturbo_log2_fc": [-1.0, 0.0],
            "perturbo_p_value": [0.01, 0.9],
        }
    )
    monkeypatch.setattr("evaluate_controls.plot_volcano", lambda *args, **kwargs: None)
    monkeypatch.setattr("evaluate_controls.perform_binary_evaluation", lambda *args, **kwargs: True)
    monkeypatch.setattr("evaluate_controls.savefig", lambda *args, **kwargs: None)

    run_evaluation_controls(mdata, outdir=tmp_path)

    assert not (tmp_path / "controls_evaluation_skipped.txt").exists()


def test_untested_pairs_do_not_unbalance_the_matched_classes(tmp_path, monkeypatch):
    import matplotlib

    matplotlib.use("Agg", force=True)

    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(
        X=np.ones((1, 2)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["ENSG1", "ENSG2"]),
    )
    guide_ids = ["t1", "t2", "c1", "c2", "c3"]
    guide_var = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "intended_target_name": [
                "GENE1",
                "GENE2",
                "non-targeting|1",
                "non-targeting|1",
                "non-targeting|2",
            ],
            "targeting": [True, True, False, False, False],
        },
        index=guide_ids,
    )
    guide = ad.AnnData(X=np.ones((1, 5)), obs=obs.copy(), var=guide_var)
    mdata = mu.MuData({"gene": gene, "guide": guide})
    # t2/ENSG2 is a direct target the conditional randomization test did not
    # test, so it carries no p-value; c3/ENSG1 is an untested control.
    mdata.uns["global_analysis_per_guide_results"] = pd.DataFrame(
        {
            "guide_id": ["t1", "t2", "c1", "c2", "c3"],
            "gene_id": ["ENSG1", "ENSG2", "ENSG1", "ENSG2", "ENSG1"],
            "gene_name": ["GENE1", "GENE2", "GENE1", "GENE2", "GENE1"],
            "perturbo_log2_fc": [-1.0, -1.0, 0.0, 0.0, 0.0],
            "perturbo_p_value": [0.01, np.nan, 0.4, 0.5, np.nan],
        }
    )
    monkeypatch.setattr("evaluate_controls.plot_volcano", lambda *args, **kwargs: None)
    monkeypatch.setattr("evaluate_controls.savefig", lambda *args, **kwargs: None)

    run_evaluation_controls(mdata, outdir=tmp_path)

    summary = (tmp_path / "controls_evaluation_summary.txt").read_text()
    assert "direct_target_rows=1" in summary
    assert "non_targeting_rows=1" in summary
    # The classes are matched after the untested pairs are removed, so nothing
    # is left for the curve builder to drop.
    assert "unscorable_rows_dropped=0" in summary
    assert "positive_rate=0.500000" in summary


def test_control_evaluation_skips_cleanly_without_global_results(tmp_path):
    class LocalOnlyResult:
        uns = {"local_analysis_per_guide_results": pd.DataFrame()}

    run_evaluation_controls(LocalOnlyResult(), outdir=tmp_path)

    marker = tmp_path / "controls_evaluation_skipped.txt"
    assert marker.exists()
    assert "global PerTurbo guide results are not present" in marker.read_text()
