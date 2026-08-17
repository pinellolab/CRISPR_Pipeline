import pathlib
import sys

import numpy as np
from scipy import sparse


BIN_DIR = pathlib.Path(__file__).resolve().parents[1] / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from qc_metrics_json import _guide_assignment_summary


class GuideMod:
    def __init__(self, assignment=None):
        self.layers = {}
        if assignment is not None:
            self.layers["guide_assignment"] = assignment


def test_guide_assignment_summary_uses_assignment_layer_for_means():
    assignment = sparse.csr_matrix([[1, 0, 1], [0, 1, 0]])
    summary = _guide_assignment_summary(GuideMod(assignment))

    assert summary["total_sgrna_assignment_values"] == 3.0
    assert summary["median_sgrna_assignment_per_guide"] == 1.0
    assert summary["mean_guides_per_cell"] == 1.5
    assert summary["mean_cells_per_guide"] == 1.0


def test_guide_assignment_summary_returns_none_without_assignment_layer():
    summary = _guide_assignment_summary(GuideMod())

    assert summary["mean_guides_per_cell"] is None
    assert summary["mean_cells_per_guide"] is None


def test_assignment_means_are_not_raw_umi_means():
    assignment = np.array([[1, 0], [0, 1]])
    raw_umis = np.array([[100, 0], [0, 200]])

    summary = _guide_assignment_summary(GuideMod(assignment))

    assert summary["mean_guides_per_cell"] == 1.0
    assert summary["mean_guides_per_cell"] != raw_umis.sum(axis=1).mean()
