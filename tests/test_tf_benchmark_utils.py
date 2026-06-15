import sys
from pathlib import Path

import pandas as pd
import pytest


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN_DIR))

from tf_benchmark_utils import normalize_tf_benchmark_results


def test_normalizes_merged_perturbo_schema():
    results = pd.DataFrame(
        {
            "gene_id": ["ENSG_TARGET"],
            "intended_target_name": ["ENSG_TF"],
            "perturbo_log2_fc": [0.5],
            "perturbo_p_value": [0.001],
        }
    )

    normalized = normalize_tf_benchmark_results(results)

    assert normalized.loc[0, "log2_fc"] == 0.5
    assert normalized.loc[0, "p_value"] == 0.001


def test_preserves_generic_schema_and_fills_missing_values_from_alias():
    results = pd.DataFrame(
        {
            "p_value": [0.01, None],
            "perturbo_p_value": [0.02, 0.03],
        }
    )

    normalized = normalize_tf_benchmark_results(results)

    assert normalized["p_value"].tolist() == [0.01, 0.03]


def test_rejects_results_without_supported_pvalue_column():
    with pytest.raises(ValueError, match="p_value, perturbo_p_value"):
        normalize_tf_benchmark_results(pd.DataFrame({"gene_id": ["ENSG1"]}))
