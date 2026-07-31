"""Memory-safe helpers for dashboard matrix summaries."""

from __future__ import annotations

import numpy as np
from scipy import sparse


def guide_frequencies(guide_assignment_matrix) -> np.ndarray:
    """Return per-guide assignment totals for dense or sparse 2-D matrices."""
    if sparse.issparse(guide_assignment_matrix):
        if guide_assignment_matrix.ndim != 2:
            raise ValueError("guide_assignment must be a two-dimensional matrix")
        return np.asarray(guide_assignment_matrix.sum(axis=0)).ravel()

    values = np.asarray(guide_assignment_matrix)
    if values.ndim != 2:
        raise ValueError("guide_assignment must be a two-dimensional matrix")
    return np.asarray(values.sum(axis=0)).ravel()
