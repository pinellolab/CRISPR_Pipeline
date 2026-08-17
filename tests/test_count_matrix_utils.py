import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from count_matrix_utils import (
    normalize_sparse_index_dtypes,
    smallest_count_dtype,
    to_sparse_counts,
)


def test_smallest_count_dtype_prefers_uint16():
    assert smallest_count_dtype(0) == np.uint16
    assert smallest_count_dtype(65535) == np.uint16


def test_smallest_count_dtype_widens_past_uint16():
    assert smallest_count_dtype(65536) == np.uint32
    assert smallest_count_dtype(5_000_000_000) == np.uint64


def test_smallest_count_dtype_uses_signed_when_negative():
    assert smallest_count_dtype(10, min_value=-5) == np.int16
    assert smallest_count_dtype(10, min_value=-40000) == np.int32


def test_dense_counts_become_csr_uint16():
    dense = np.array([[0, 3, 0], [7, 0, 0]], dtype=np.float64)
    observed = to_sparse_counts(dense)
    assert sparse.isspmatrix_csr(observed)
    assert observed.dtype == np.uint16
    np.testing.assert_array_equal(observed.toarray(), dense)


def test_no_silent_wraparound_above_uint16_max():
    """A blind .astype(np.uint16) turns 70000 into 4464; this must not."""
    dense = np.array([[65535, 70000]], dtype=np.int64)
    observed = to_sparse_counts(dense)
    assert observed.dtype == np.uint32
    np.testing.assert_array_equal(observed.toarray(), dense)


def test_non_integral_data_is_preserved_as_float_not_truncated():
    dense = np.array([[0.0, 1.5], [2.25, 0.0]])
    observed = to_sparse_counts(dense)
    assert observed.dtype == np.float32
    np.testing.assert_allclose(observed.toarray(), dense)


def test_negative_integral_data_uses_signed_dtype():
    dense = np.array([[-3, 0, 5]], dtype=np.int64)
    observed = to_sparse_counts(dense)
    assert observed.dtype == np.int16
    np.testing.assert_array_equal(observed.toarray(), dense)


def test_accepts_pandas_dataframe():
    df = pd.DataFrame({"g1": [0, 1], "g2": [1, 0]})
    observed = to_sparse_counts(df)
    assert sparse.isspmatrix_csr(observed)
    assert observed.dtype == np.uint16
    np.testing.assert_array_equal(observed.toarray(), df.to_numpy())


def test_accepts_existing_sparse_and_preserves_values():
    original = sparse.coo_matrix(np.array([[0, 2], [9, 0]], dtype=np.float64))
    observed = to_sparse_counts(original)
    assert sparse.isspmatrix_csr(observed)
    assert observed.dtype == np.uint16
    np.testing.assert_array_equal(observed.toarray(), original.toarray())


def test_empty_and_all_zero_matrices():
    assert to_sparse_counts(np.zeros((3, 3))).dtype == np.uint16
    empty = to_sparse_counts(sparse.csr_matrix((2, 2)))
    assert empty.nnz == 0
    assert empty.dtype == np.uint16


def test_boolean_assignment_matrix_becomes_uint16():
    dense = np.array([[True, False], [False, True]])
    observed = to_sparse_counts(dense)
    assert observed.dtype == np.uint16
    np.testing.assert_array_equal(observed.toarray(), dense.astype(np.uint16))


def test_normalize_sparse_index_dtypes_repairs_scipy_compressed_operations():
    """Reproduce the large concat failure seen in PreprocessAnnData.

    SciPy rejects a CSR matrix whose indices and indptr use different integer
    widths. The repair must leave compact count data untouched while making
    eliminate_zeros (used by Scanpy QC) valid again.
    """
    matrix = sparse.csr_matrix(
        np.array([[0, 3, 0], [7, 0, 0]], dtype=np.uint16)
    )
    matrix.indices = matrix.indices.astype(np.int32)
    matrix.indptr = matrix.indptr.astype(np.int64)

    with pytest.raises(ValueError, match="Output dtype not compatible"):
        matrix.eliminate_zeros()

    observed = normalize_sparse_index_dtypes(matrix)
    assert observed is matrix
    assert observed.dtype == np.uint16
    assert observed.indices.dtype == observed.indptr.dtype == np.dtype(np.int32)
    observed.eliminate_zeros()
    np.testing.assert_array_equal(
        observed.toarray(), np.array([[0, 3, 0], [7, 0, 0]], dtype=np.uint16)
    )


def test_normalize_sparse_index_dtypes_is_noop_when_already_compatible():
    matrix = sparse.csr_matrix(np.eye(3, dtype=np.uint16))
    original_indices = matrix.indices
    original_indptr = matrix.indptr

    observed = normalize_sparse_index_dtypes(matrix)

    assert observed is matrix
    assert observed.indices is original_indices
    assert observed.indptr is original_indptr


@pytest.mark.parametrize("axis", [0, 1])
def test_sums_do_not_overflow_after_narrowing(axis):
    """uint16 storage must not corrupt per-cell/per-feature totals, which
    routinely exceed 65535 even when no individual count does."""
    n = 1000
    dense = np.full((n, 3), 1000, dtype=np.int64)  # totals = 1e6 >> 65535
    stored = to_sparse_counts(dense)
    assert stored.dtype == np.uint16

    observed = np.asarray(stored.sum(axis=axis)).ravel()
    expected = dense.sum(axis=axis)
    np.testing.assert_array_equal(observed, expected)


def test_sum_with_explicit_int64_dtype_matches():
    dense = np.full((500, 2), 900, dtype=np.int64)
    stored = to_sparse_counts(dense)
    observed = np.asarray(stored.sum(axis=0, dtype=np.int64)).ravel()
    np.testing.assert_array_equal(observed, dense.sum(axis=0))


def test_fractional_mapper_counts_are_not_truncated():
    """kb-python's nac workflow (cc-perturb-seq) EM-distributes ambiguous
    reads, so mature/nascent/ambiguous layers can hold fractional counts.
    Narrowing those to an integer dtype would silently floor real signal --
    e.g. 0.5 -> 0. They must stay float."""
    ambiguous = np.array([[0.0, 0.5], [1.25, 0.0]])
    observed = to_sparse_counts(ambiguous)
    assert observed.dtype == np.float32
    np.testing.assert_allclose(observed.toarray(), ambiguous)
    assert observed.toarray()[0, 1] == pytest.approx(0.5)


def test_summed_nac_layers_stay_float_when_fractional():
    """mature + nascent + ambiguous is what lands in .X for the nac
    workflow; if any component is fractional the sum must remain float."""
    mature = np.array([[1.0, 2.0]])
    nascent = np.array([[0.0, 1.0]])
    ambiguous = np.array([[0.5, 0.25]])
    combined = mature + nascent + ambiguous

    observed = to_sparse_counts(combined)
    assert observed.dtype == np.float32
    np.testing.assert_allclose(observed.toarray(), combined)


def test_float_valued_but_whole_counts_still_narrow_to_int():
    """The common case: mapper emits float64 that happens to hold whole
    numbers (e.g. after rounding). Those are safe to store as uint16."""
    whole_floats = np.array([[0.0, 3.0], [7.0, 0.0]], dtype=np.float64)
    observed = to_sparse_counts(whole_floats)
    assert observed.dtype == np.uint16
    np.testing.assert_array_equal(observed.toarray(), whole_floats)


def test_nan_and_inf_data_is_not_coerced_to_int():
    """NaN/inf can't be represented as an integer; casting would produce
    garbage rather than failing, so such data must stay float."""
    with_nan = np.array([[0.0, np.nan], [1.0, 0.0]])
    observed = to_sparse_counts(with_nan)
    assert observed.dtype == np.float32
    assert np.isnan(observed.toarray()[0, 1])


def test_anndata_roundtrip_preserves_values_and_dtype(tmp_path):
    """Values must survive an actual h5ad write/read at the narrowed dtype."""
    ad = pytest.importorskip("anndata")
    dense = np.array([[0, 5, 0], [70000, 0, 2]], dtype=np.int64)
    adata = ad.AnnData(
        X=to_sparse_counts(dense),
        obs=pd.DataFrame(index=["c1", "c2"]),
        var=pd.DataFrame(index=["g1", "g2", "g3"]),
    )
    adata.layers["guide_assignment"] = to_sparse_counts((dense > 0).astype(int))

    path = tmp_path / "roundtrip.h5ad"
    adata.write_h5ad(path)
    reloaded = ad.read_h5ad(path)

    assert sparse.issparse(reloaded.X)
    np.testing.assert_array_equal(reloaded.X.toarray(), dense)
    np.testing.assert_array_equal(
        reloaded.layers["guide_assignment"].toarray(), (dense > 0).astype(np.uint16)
    )
    assert reloaded.layers["guide_assignment"].dtype == np.uint16
