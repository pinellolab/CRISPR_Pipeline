import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pytest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import add_perturbo_results_to_mudata
from analysis_output_formatting import make_h5mu_safe_dataframe
from result_table_io import read_result_table, write_result_table


def _base_mudata():
    obs = pd.DataFrame(index=["cell1"])
    gene = ad.AnnData(
        X=np.ones((1, 1)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["GENE1"]),
    )
    guide = ad.AnnData(
        X=np.ones((1, 1)),
        obs=obs.copy(),
        var=pd.DataFrame(index=["guide1"]),
    )
    return mu.MuData({"gene": gene, "guide": guide})


def test_parquet_result_roundtrip(tmp_path):
    pytest.importorskip("pyarrow")
    expected = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "p_value": [0.01, 0.2],
            "perturbo_q_value": [0.02, 0.2],
        }
    )
    output = tmp_path / "results.parquet"

    write_result_table(expected, output)
    observed = read_result_table(output)

    pd.testing.assert_frame_equal(observed, expected, check_dtype=False)


def test_add_results_normalizes_mixed_chromosome_values_for_h5mu(tmp_path):
    pytest.importorskip("pyarrow")
    guide_results = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE1"],
            "guide_id": ["guide1", "guide2"],
            "log2_fc": [0.1, 0.2],
            "p_value": [0.5, 0.4],
        }
    )
    element_results = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE1", "GENE1"],
            "intended_target_name": ["target1", "target2", "target3"],
            "intended_target_chr": ["chr8", 8, None],
            "intended_target_start": [1, 2, 3],
            "intended_target_end": [2, 3, 4],
            "log2_fc": [0.1, 0.2, 0.3],
            "p_value": [0.5, 0.4, 0.3],
        }
    )
    guide_path = tmp_path / "guide.parquet"
    element_path = tmp_path / "element.parquet"
    base_path = tmp_path / "base.h5mu"
    output_path = tmp_path / "inference.h5mu"
    write_result_table(guide_results, guide_path)
    write_result_table(element_results, element_path)
    _base_mudata().write(base_path)

    add_perturbo_results_to_mudata.add_perturbo_results_to_mudata(
        guide_path,
        element_path,
        base_path,
        output_path,
    )

    observed = mu.read_h5mu(output_path)
    stored = pd.DataFrame(observed.uns["per_element_results"])
    assert stored["intended_target_chr"].tolist() == ["chr8", "8", ""]


def test_h5mu_safe_dataframe_stringifies_mixed_objects():
    frame = pd.DataFrame({"mixed": ["chr8", 8, None]})
    observed = make_h5mu_safe_dataframe(frame)
    assert observed["mixed"].tolist() == ["chr8", "8", ""]
