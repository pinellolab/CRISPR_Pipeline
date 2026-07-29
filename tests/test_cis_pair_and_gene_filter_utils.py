import importlib.util
from pathlib import Path
import sys


BIN_DIR = Path(__file__).resolve().parents[1] / "bin"


def load_module(name):
    sys.path.insert(0, str(BIN_DIR))
    spec = importlib.util.spec_from_file_location(name, BIN_DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_chromosome_namespace_normalization():
    pairs = load_module("create_pairs_to_test")
    assert pairs.normalize_chromosome("chr8") == "8"
    assert pairs.normalize_chromosome("8") == "8"
    assert pairs.normalize_chromosome("chrX") == "X"
    assert pairs.normalize_chromosome("x") == "X"
    assert pairs.normalize_chromosome(None) is None


def test_fractional_gene_support_preserves_historical_boundary():
    concat = load_module("mudata_concat")
    assert concat.resolve_min_cells(126_154, 0.05) == 6_308
    assert concat.resolve_min_cells(100, 0) == 1


def test_absolute_gene_support_is_rejected():
    concat = load_module("mudata_concat")
    try:
        concat.resolve_min_cells(126_154, 15)
    except ValueError as error:
        assert "fraction" in str(error)
    else:
        raise AssertionError("Absolute gene-support thresholds must be rejected")


def test_run_param_preflight_rejects_stale_absolute_threshold():
    validator = load_module("validate_fractional_qc_params")
    value, tapseq_mode = validator.validate_fractional_qc_params(
        {"QC_min_cells_per_gene": 0.000018, "TAPSEQ_QC_MODE": False}
    )
    assert value == 0.000018
    assert tapseq_mode is False

    try:
        validator.validate_fractional_qc_params({"QC_min_cells_per_gene": 15})
    except ValueError as error:
        assert "fraction" in str(error)
    else:
        raise AssertionError("Stale absolute thresholds must fail preflight")
