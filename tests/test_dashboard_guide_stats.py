#!/usr/bin/env python
"""Regression tests for the dashboard guide-assignment summary metrics.

Run with pytest, or standalone with the pipeline's own dependencies:

    python tests/test_dashboard_guide_stats.py
"""

import os
import sys

import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'bin')
)

from create_dashboard_df import guide_assignment_counts
from create_dashboard_df_HASHING import (
    guide_assignment_counts as guide_assignment_counts_hashing,
)

# Three cells x four guides. The raw guide UMI counts carry background signal on
# almost every guide, while the assignment matrix keeps one guide for cell 0, two
# for cell 1 and none for cell 2.
RAW_GUIDE_UMIS = np.array([
    [50, 3, 1, 0],
    [40, 0, 30, 2],
    [1, 2, 0, 1],
], dtype=float)

ASSIGNMENT = np.array([
    [1, 0, 0, 0],
    [1, 0, 1, 0],
    [0, 0, 0, 0],
], dtype=float)

IMPLEMENTATIONS = (guide_assignment_counts, guide_assignment_counts_hashing)


def make_guide_mod(assignment, raw_counts=RAW_GUIDE_UMIS, sparse=False):
    """A guide modality shaped like the one add_guide_assignment.py writes."""
    guide = ad.AnnData(X=csr_matrix(raw_counts) if sparse else raw_counts.copy())
    guide.layers['guide_assignment'] = (
        csr_matrix(assignment) if sparse else assignment.copy()
    )
    return guide


def test_guides_per_cell_counts_assigned_guides_not_umis():
    for counts in IMPLEMENTATIONS:
        for sparse in (False, True):
            guides_per_cell = counts(make_guide_mod(ASSIGNMENT, sparse=sparse), axis=1)
            assert list(guides_per_cell) == [1, 2, 0]
            assert np.mean(guides_per_cell) == 1.0
            # Summing .X instead would report 42.33 guides per cell.
            assert np.mean(guides_per_cell) != np.mean(np.sum(RAW_GUIDE_UMIS, axis=1))


def test_cells_per_guide_counts_assigned_cells_not_umis():
    for counts in IMPLEMENTATIONS:
        for sparse in (False, True):
            cells_per_guide = counts(make_guide_mod(ASSIGNMENT, sparse=sparse), axis=0)
            assert list(cells_per_guide) == [2, 0, 1, 0]
            assert np.mean(cells_per_guide) == 0.75
            assert np.mean(cells_per_guide) != np.mean(np.sum(RAW_GUIDE_UMIS, axis=0))


def test_summary_follows_the_assignment_method():
    """Two assignment callers over identical raw counts must summarise differently."""
    permissive = make_guide_mod((RAW_GUIDE_UMIS > 0).astype(float))
    strict = make_guide_mod((RAW_GUIDE_UMIS >= 30).astype(float))
    for counts in IMPLEMENTATIONS:
        assert np.mean(counts(permissive, axis=1)) != np.mean(counts(strict, axis=1))
        assert np.mean(counts(permissive, axis=0)) != np.mean(counts(strict, axis=0))


def test_non_binary_assignment_values_count_once_per_guide():
    """An assignment layer holding scores rather than 0/1 still counts guides."""
    for counts in IMPLEMENTATIONS:
        assert list(counts(make_guide_mod(ASSIGNMENT * 7.0), axis=1)) == [1, 2, 0]
        assert list(
            counts(make_guide_mod(ASSIGNMENT * 7.0, sparse=True), axis=0)
        ) == [2, 0, 1, 0]


def test_empty_assignment_layer():
    empty = np.zeros((3, 4))
    for counts in IMPLEMENTATIONS:
        assert list(counts(make_guide_mod(empty), axis=1)) == [0, 0, 0]
        assert list(counts(make_guide_mod(empty, sparse=True), axis=0)) == [0, 0, 0, 0]


if __name__ == '__main__':
    for name, test in list(globals().items()):
        if name.startswith('test_') and callable(test):
            test()
            print('PASS', name)
