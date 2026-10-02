"""Storage rounding must not change significance probabilities or raw inputs."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from compact_result_dtypes import compact_result_floats


def test_compact_effects_preserves_probabilities_and_input(tmp_path):
    values = np.array([1.234567891, -2.123456789, 0, np.nan, np.inf])
    probabilities = np.array([1e-300, np.nextafter(0.0, 1.0), 0.05, np.nan, 1.0])
    original = pd.DataFrame({
        "log2_fc": values, "perturbo_fc_se": np.abs(values),
        "perturbo_negLog10p": [300., 323.306215343, 1.301029996, np.nan, 0.],
        "p_value": probabilities, "perturbo_q_value": probabilities,
        "perturbo_posterior_prob": probabilities,
        "unknown_metric": values, "nPerturbedCells": [1, 2, 3, 4, 5],
        "guide_id": pd.Categorical(["a", "b", "a", "b", "a"]),
    })
    before = original.copy(deep=True)
    packed = compact_result_floats(original)
    pd.testing.assert_frame_equal(original, before)
    for name in ("log2_fc", "perturbo_fc_se", "perturbo_negLog10p"):
        assert packed[name].dtype == np.float32
        np.testing.assert_allclose(packed[name], original[name], rtol=6e-8, atol=0,
                                   equal_nan=True)
    for name in original.columns.difference(["log2_fc", "perturbo_fc_se", "perturbo_negLog10p"]):
        pd.testing.assert_series_equal(packed[name], original[name])
    path = tmp_path / "compact.parquet"
    packed.to_parquet(path, index=False)
    pd.testing.assert_frame_equal(pd.read_parquet(path), packed)
    assert packed.memory_usage(deep=True).sum() < original.memory_usage(deep=True).sum()


@pytest.mark.parametrize("value", [1e100, -1e100, 1e-100, -1e-100])
def test_unrepresentable_effect_fails_without_corrupting_raw(value):
    source = pd.DataFrame({"log2_fc": [value]})
    with pytest.raises(ValueError, match="overflow or underflow"):
        compact_result_floats(source)
    assert source.loc[0, "log2_fc"] == value
