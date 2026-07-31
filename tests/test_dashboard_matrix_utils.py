import pathlib
import sys

import numpy as np
from scipy import sparse


BIN_DIR = pathlib.Path(__file__).resolve().parents[1] / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from dashboard_matrix_utils import guide_frequencies


def test_guide_frequencies_accepts_dense_matrix():
    matrix = np.array([[1, 0, 2], [0, 3, 1]])

    np.testing.assert_array_equal(guide_frequencies(matrix), [1, 3, 3])


def test_guide_frequencies_accepts_sparse_matrix():
    matrix = sparse.csr_matrix([[1, 0, 2], [0, 3, 1]])

    np.testing.assert_array_equal(guide_frequencies(matrix), [1, 3, 3])


def test_guide_frequencies_rejects_non_matrix():
    with np.testing.assert_raises_regex(ValueError, "two-dimensional"):
        guide_frequencies(np.array([1, 2, 3]))
