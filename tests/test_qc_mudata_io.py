"""The slim readers must agree with read_h5mu, column for column.

They exist to stop the QC scripts loading the full result tables out of uns --
6.35 GB of a 6.41 GB file on the TAP-seq chr8 screen -- so the thing to pin is
that reading less returns the same values, not merely that it is faster.
"""

import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import qc_mudata_io as qc

KEY = "global_analysis_per_guide_results"
N = 12


def _results_frame():
    rng = np.random.default_rng(3)
    return pd.DataFrame(
        {
            "gene_id": [f"ENSG{i:04d}" for i in range(N)],
            "guide_id": [f"sg{i:03d}" for i in range(N)],
            # categorical: how mergedResults stores low-cardinality text
            "intended_target_name": pd.Categorical(
                ["elemA", "elemB"] * (N // 2)
            ),
            "gene_name": pd.Categorical([f"SYM{i%3}" for i in range(N)]),
            "targeting": [True, False] * (N // 2),
            "log2_fc": rng.normal(size=N),
            "p_value": np.concatenate([rng.random(N - 2), [0.0, 1e-320]]),
            "q_value": rng.random(N),
            # nullable integer, as intended_target_start arrives
            "nPerturbedCells": pd.array(
                [10, 20, None] * (N // 3), dtype="Int64"
            ),
            # a heavy column the QC scripts never read
            "unused_metric": rng.normal(size=N),
        }
    )


def _mudata(tmp_path):
    obs = pd.DataFrame(index=[f"cell{i}" for i in range(6)])
    gene = ad.AnnData(
        X=sparse.csr_matrix(np.arange(12, dtype=np.float32).reshape(6, 2)),
        obs=obs.copy(),
        var=pd.DataFrame({"symbol": ["S1", "S2"]}, index=["G1", "G2"]),
    )
    guide = ad.AnnData(
        X=sparse.csr_matrix(np.eye(6, 3, dtype=np.float32)),
        obs=obs.copy(),
        var=pd.DataFrame(
            {"guide_id": ["a", "b", "c"], "targeting": [True, True, False]},
            index=["a", "b", "c"],
        ),
    )
    guide.layers["guide_assignment"] = guide.X.copy()
    mdata = mu.MuData({"gene": gene, "guide": guide})
    mdata.uns[KEY] = _results_frame()
    mdata.uns["local_analysis_per_guide_results"] = _results_frame().head(3)
    path = tmp_path / "inference_mudata.h5mu"
    mdata.write(path)
    return path


def test_result_keys_lists_tables_without_reading(tmp_path):
    path = _mudata(tmp_path)
    assert qc.result_keys(path) == [
        "global_analysis_per_guide_results",
        "local_analysis_per_guide_results",
    ]


def test_result_table_columns_without_reading(tmp_path):
    path = _mudata(tmp_path)
    got = set(qc.result_table_columns(path, KEY))
    assert {"gene_id", "guide_id", "p_value", "log2_fc"} <= got


@pytest.mark.parametrize(
    "column",
    ["gene_id", "guide_id", "intended_target_name", "gene_name", "targeting",
     "log2_fc", "p_value", "q_value", "nPerturbedCells"],
)
def test_read_result_columns_matches_read_h5mu(tmp_path, column):
    path = _mudata(tmp_path)
    full = mu.read_h5mu(path).uns[KEY]
    slim = qc.read_result_columns(path, KEY, [column])
    assert list(slim.columns) == [column]
    expected = pd.Series(full[column]).reset_index(drop=True)
    actual = pd.Series(slim[column]).reset_index(drop=True)
    # Compare as objects: the slim reader returns plain arrays where read_h5mu
    # may hand back pandas extension dtypes carrying the same values.
    e = [None if pd.isna(v) else v for v in expected.tolist()]
    a = [None if pd.isna(v) else v for v in actual.tolist()]
    assert a == e, column


def test_read_result_columns_skips_absent_and_unrequested(tmp_path):
    path = _mudata(tmp_path)
    slim = qc.read_result_columns(path, KEY, ["gene_id", "not_a_column"])
    assert list(slim.columns) == ["gene_id"]
    assert "unused_metric" not in slim.columns


def test_read_mudata_without_uns_keeps_modalities_drops_tables(tmp_path):
    path = _mudata(tmp_path)
    full = mu.read_h5mu(path)
    slim = qc.read_mudata_without_uns(path)
    assert set(slim.mod) == set(full.mod)
    assert not any(k in slim.uns for k in qc.RESULTS_KEY_CANDIDATES)
    for name in full.mod:
        assert slim[name].shape == full[name].shape
        pd.testing.assert_frame_equal(slim[name].var, full[name].var)
        np.testing.assert_allclose(
            slim[name].X.toarray() if sparse.issparse(slim[name].X) else slim[name].X,
            full[name].X.toarray() if sparse.issparse(full[name].X) else full[name].X,
        )
    assert "guide_assignment" in slim["guide"].layers


def test_resolve_result_key_prefers_requested_then_candidates(tmp_path):
    path = _mudata(tmp_path)
    assert qc.resolve_result_key(path, KEY, qc.RESULTS_KEY_CANDIDATES) == KEY
    assert qc.resolve_result_key(path, "auto", qc.RESULTS_KEY_CANDIDATES) == KEY
    # A requested-but-absent key falls back rather than failing.
    assert qc.resolve_result_key(path, "nope", qc.RESULTS_KEY_CANDIDATES) == KEY
    assert qc.resolve_result_key(path, "auto", ("only_missing",)) is None
