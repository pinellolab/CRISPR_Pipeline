import importlib.util
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "bin" / "run_gmm_demux_with_qc.py"
SPEC = importlib.util.spec_from_file_location("run_gmm_demux_with_qc", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_rejects_zero_versus_nonzero_fit():
    counts = pd.DataFrame(
        {"HTO1": [0, 1, 2, 10], "HTO2": [0, 0, 8, 1]},
        index=["a", "b", "c", "d"],
    )
    assignments = pd.DataFrame(
        {
            "Cluster_id": [0, 1, 3, 3],
            "Confidence": [0.9, 0.4, 0.6, 0.7],
            "hto_type": ["negative", "HTO1", "HTO1-HTO2", "HTO1-HTO2"],
        },
        index=counts.index,
    )

    result = MODULE.evaluate_fit(counts, assignments, reject_nonzero_positive=True)

    assert not result["accepted"]
    assert result["per_hto"]["HTO1"]["exact_nonzero_equivalence"]
    assert result["per_hto"]["HTO1"]["rejection_reason"]


def test_accepts_positive_component_with_nonzero_background():
    counts = pd.DataFrame(
        {"HTO1": [0, 1, 2, 6, 10], "HTO2": [0, 2, 0, 1, 8]},
        index=["a", "b", "c", "d", "e"],
    )
    assignments = pd.DataFrame(
        {
            "Cluster_id": [0, 0, 0, 1, 2],
            "Confidence": [0.9, 0.8, 0.8, 0.9, 0.9],
            "hto_type": ["negative", "negative", "negative", "HTO1", "HTO2"],
        },
        index=counts.index,
    )

    result = MODULE.evaluate_fit(counts, assignments, reject_nonzero_positive=True)

    assert result["accepted"]
    assert result["per_hto"]["HTO1"]["minimum_positive_count"] == 6
    assert not result["per_hto"]["HTO1"]["near_nonzero_equivalence"]


def test_gate_can_be_disabled_explicitly():
    counts = pd.DataFrame({"HTO1": [0, 1]}, index=["a", "b"])
    assignments = pd.DataFrame(
        {"Cluster_id": [0, 1], "Confidence": [0.9, 0.9], "hto_type": ["negative", "HTO1"]},
        index=counts.index,
    )

    result = MODULE.evaluate_fit(counts, assignments, reject_nonzero_positive=False)

    assert result["accepted"]
    assert result["per_hto"]["HTO1"]["exact_nonzero_equivalence"]
