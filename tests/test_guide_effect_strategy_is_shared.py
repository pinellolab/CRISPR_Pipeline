"""The pipeline must keep running shared guide effects across the PerTurbo upgrade.

PerTurbo v1 called the continuous per-guide efficacy ``scaled``. v2 renamed it
``relative`` and added ``shared``, which fits none - but through 2.0.0rc11 the
Python API silently mapped ``scaled`` onto ``shared``. So every run of this
pipeline to date has used shared guide effects while asking for ``scaled``. From
rc12 that name resolves to ``relative`` instead, which would change the model
under a container bump, with no diff in this repository to show for it.

These tests pin the strategy to shared at the only place it is decided.
"""

import argparse
import pathlib
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

pytest.importorskip("perturbo")
pytest.importorskip("scvi")

import perturbo_inference as pi  # noqa: E402


def test_every_accepted_efficiency_mode_resolves_to_shared() -> None:
    assert pi.GUIDE_EFFECT_STRATEGY == "shared"
    for value in ("scaled", "shared", None):
        assert pi.resolve_efficiency_mode(value) == "shared"


def test_a_mode_that_would_fit_efficacy_is_refused() -> None:
    with pytest.raises(ValueError, match="not supported"):
        pi.resolve_efficiency_mode("relative")


def test_the_deprecated_boolean_flag_parses_false_as_false() -> None:
    """argparse's type=bool maps every non-empty string to True, so
    ``--fit_guide_efficacy False`` used to arrive as True."""
    assert pi._parse_bool("False") is False
    assert pi._parse_bool("false") is False
    assert pi._parse_bool("0") is False
    assert pi._parse_bool("True") is True
    with pytest.raises(argparse.ArgumentTypeError):
        pi._parse_bool("maybe")


def test_the_chunked_driver_does_not_forward_the_deprecated_flag() -> None:
    import perturbo_inference_chunked as pic

    command = pic.build_perturbo_command(
        perturbo_script=BIN_DIR / "perturbo_inference.py",
        mdata_input_fp="in.h5mu",
        results_tsv_fp="out.tsv",
    )
    assert "--fit_guide_efficacy" not in command
    assert command[command.index("--efficiency_mode") + 1] == "shared"
