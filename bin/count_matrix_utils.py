#!/usr/bin/env python
"""Helpers for storing count matrices compactly in AnnData/MuData.

Single-cell count matrices (RNA UMIs, guide UMIs, guide assignments) are
overwhelmingly zeros and hold small non-negative integers, but are frequently
built as dense float64 arrays -- costing ~4x the memory/disk of a uint16 CSR
matrix before sparsity is even considered.

The dtype is always chosen from the data's actual maximum, never assumed:
``np.uint16`` covers per-cell/per-gene counts in practice, but a blind
``.astype(np.uint16)`` silently *wraps* anything above 65535 (70000 -> 4464)
rather than raising, which would corrupt counts for deeply sequenced cells or
highly expressed genes. ``smallest_count_dtype`` widens to uint32/uint64
instead.

Summation is safe in both directions: numpy and scipy both upcast uint16 to
uint64 for ``.sum()``, so per-cell/per-gene totals do not overflow. Callers
that narrow a sum back down should pass an explicit ``dtype``.
"""

import numpy as np
from scipy import sparse

# Narrowest first. uint16 (max 65535) covers the overwhelming majority of
# per-cell/per-feature UMI counts; wider types are fallbacks, not defaults.
UNSIGNED_COUNT_DTYPES = (np.uint16, np.uint32, np.uint64)
SIGNED_COUNT_DTYPES = (np.int16, np.int32, np.int64)


def smallest_count_dtype(max_value, min_value=0):
    """Return the narrowest integer dtype that holds [min_value, max_value].

    Uses unsigned types when the data is non-negative (the normal case for
    counts), signed types otherwise.
    """
    candidates = UNSIGNED_COUNT_DTYPES if min_value >= 0 else SIGNED_COUNT_DTYPES
    for dtype in candidates:
        info = np.iinfo(dtype)
        if min_value >= info.min and max_value <= info.max:
            return np.dtype(dtype)
    return np.dtype(candidates[-1])


def _matrix_extremes(matrix):
    """Return (min, max) over stored values, treating an empty matrix as (0, 0)."""
    if sparse.issparse(matrix):
        if matrix.nnz == 0:
            return 0, 0
        data = matrix.data
        # A sparse matrix with implicit zeros can't have a positive minimum.
        return min(0, data.min()), data.max()
    values = np.asarray(matrix)
    if values.size == 0:
        return 0, 0
    return values.min(), values.max()


def _is_integral(matrix):
    """True if every stored value is a whole, finite number."""
    data = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
    if data.size == 0:
        return True
    if np.issubdtype(data.dtype, np.integer) or data.dtype == np.bool_:
        return True
    if not np.issubdtype(data.dtype, np.floating):
        return False
    if not np.all(np.isfinite(data)):
        return False
    return bool(np.all(np.mod(data, 1) == 0))


def to_sparse_counts(matrix, float_dtype=np.float32):
    """Return ``matrix`` as CSR with the narrowest dtype that holds it exactly.

    Integral data (raw counts, guide assignments) becomes the narrowest
    unsigned/signed integer type -- typically uint16. Non-integral data
    (already normalized or log-transformed) is preserved as ``float_dtype``
    rather than being truncated to an integer.

    Accepts dense arrays, pandas DataFrames, and any scipy sparse format.
    """
    if hasattr(matrix, "sparse") and hasattr(matrix, "to_coo"):
        # pandas DataFrame with sparse dtypes
        matrix = matrix.sparse.to_coo()
    elif hasattr(matrix, "to_numpy"):
        matrix = matrix.to_numpy()

    if _is_integral(matrix):
        low, high = _matrix_extremes(matrix)
        target_dtype = smallest_count_dtype(int(high), int(low))
    else:
        target_dtype = np.dtype(float_dtype)

    if sparse.issparse(matrix):
        return matrix.tocsr().astype(target_dtype, copy=False)
    return sparse.csr_matrix(np.asarray(matrix), dtype=target_dtype)


def describe_matrix(matrix, label):
    """One-line summary of a stored matrix, for run logs."""
    kind = "sparse" if sparse.issparse(matrix) else "dense"
    if sparse.issparse(matrix):
        total = matrix.shape[0] * matrix.shape[1]
        density = (matrix.nnz / total) if total else 0.0
        return (
            f"{label}: {kind} {matrix.shape} dtype={matrix.dtype} "
            f"nnz={matrix.nnz} density={density:.4f}"
        )
    return f"{label}: {kind} {np.shape(matrix)} dtype={np.asarray(matrix).dtype}"
