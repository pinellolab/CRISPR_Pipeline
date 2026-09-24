import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from scipy.stats import false_discovery_control


BIN = pathlib.Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
import bounded_perturbo_results as bounded


def test_bh_sidecar_matches_scipy_for_one_complete_family(tmp_path):
    p = np.array([0.05, np.nan, 0.0, 0.05, 1.0, 1e-250, 0.2,
                  np.nextafter(0.0, 1.0), np.inf, -np.inf])
    source = tmp_path / "family.parquet"
    pd.DataFrame({"p_value": p}).to_parquet(source, row_group_size=2, index=False)

    output = bounded.write_bh_sidecar(source, tmp_path / "q.parquet")
    observed = pd.read_parquet(output).sort_values("_family_row_id")
    expected = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    expected[valid] = false_discovery_control(p[valid], method="bh")

    assert observed["perturbo_q_value"].dtype == np.float64
    np.testing.assert_array_equal(observed["_family_row_id"], np.arange(len(p)))
    np.testing.assert_allclose(
        observed["perturbo_q_value"], expected, rtol=0, atol=0, equal_nan=True
    )


def test_bh_sidecar_enforces_numeric_memory_budget(tmp_path):
    source = tmp_path / "family.parquet"
    pd.DataFrame({"p_value": [0.1, 0.2]}).to_parquet(source, index=False)
    with pytest.raises(MemoryError, match="configured"):
        bounded.write_bh_sidecar(source, tmp_path / "q.parquet", max_working_bytes=1)


def test_bh_sidecar_rejects_finite_out_of_range_values(tmp_path):
    for value in [-0.1, 1.1]:
        source = tmp_path / "family.parquet"
        pd.DataFrame({"p_value": [value]}).to_parquet(source, index=False)
        with pytest.raises(ValueError, match="p-values"):
            bounded.write_bh_sidecar(source, tmp_path / "q.parquet")


def test_manifest_detects_changed_raw_artifact(tmp_path):
    source = tmp_path / "raw.parquet"
    pd.DataFrame({"p_value": [0.1, 0.2]}).to_parquet(source, index=False)
    payload = {"raw_parquets": [bounded.parquet_manifest(source)]}
    bounded.write_manifest(tmp_path / "manifest.json", payload)
    bounded.verify_manifest(tmp_path, payload)

    pd.DataFrame({"p_value": [0.1]}).to_parquet(source, index=False)
    with pytest.raises(ValueError, match="provenance verification"):
        bounded.verify_manifest(tmp_path, payload)


@pytest.mark.parametrize("p", [[], [np.nan, np.inf, -np.inf]])
def test_empty_or_entirely_nonfinite_family(tmp_path, p):
    source = tmp_path / "family.parquet"
    pd.DataFrame({"p_value": pd.Series(p, dtype="float64")}).to_parquet(source, index=False)
    output = bounded.write_bh_sidecar(source, tmp_path / "q.parquet")
    frame = pd.read_parquet(output)
    assert len(frame) == len(p)
    assert frame.perturbo_q_value.isna().all()
