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

anndata = pytest.importorskip("anndata")
mudata = pytest.importorskip("mudata")

import collapse_guides as cg


def _make_dual_guide_mudata():
    """4 cells, each assigned exactly 2 guides targeting the same element
    (the case collapse_guides is meant to collapse), 2 elements total."""
    guide_ids = ["g1a", "g1b", "g2a", "g2b"]
    guide_var = pd.DataFrame(
        {
            "guide_id": guide_ids,
            "spacer": ["AAAA", "AAAT", "CCCC", "CCCT"],
            "targeting": [True, True, True, True],
            "type": ["targeting"] * 4,
            "pam": ["NGG"] * 4,
            "gene_name": ["GENE1", "GENE1", "GENE2", "GENE2"],
            "intended_target_name": ["elem1", "elem1", "elem2", "elem2"],
            "intended_target_chr": ["chr1", "chr1", "chr2", "chr2"],
            "intended_target_start": [100, 100, 200, 200],
            "intended_target_end": [150, 150, 250, 250],
        },
        index=guide_ids,
    )

    cell_ids = ["cell1", "cell2", "cell3", "cell4"]
    # cell1, cell2 -> both guides of elem1; cell3, cell4 -> both guides of elem2
    assignment = np.array(
        [
            [1, 1, 0, 0],
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 1, 1],
        ],
        dtype=np.uint16,
    )
    guide_x = sparse.csr_matrix(assignment)

    guide_adata = anndata.AnnData(
        X=guide_x.copy(),
        obs=pd.DataFrame(index=cell_ids),
        var=guide_var,
    )
    guide_adata.layers["guide_assignment"] = guide_x.copy()

    gene_adata = anndata.AnnData(
        X=sparse.csr_matrix(np.ones((4, 2), dtype=np.uint16)),
        obs=pd.DataFrame(index=cell_ids),
        var=pd.DataFrame(index=["GENE1", "GENE2"]),
    )

    return mudata.MuData({"gene": gene_adata, "guide": guide_adata})


def test_collapse_guides_produces_sparse_narrow_dtype_matrix(tmp_path):
    mdata = _make_dual_guide_mudata()
    input_path = tmp_path / "input.h5mu"
    mdata.write(input_path)

    output_path = tmp_path / "output.h5mu"
    result = cg.collapse_guides(str(input_path), str(output_path))

    collapsed_x = result["guide"].X
    assert sparse.issparse(collapsed_x), "collapsed guide matrix must be sparse, not dense"
    assert collapsed_x.dtype == np.uint16

    # One-hot by construction: every cell has exactly one collapsed-guide hit.
    row_sums = np.asarray(collapsed_x.sum(axis=1)).ravel()
    np.testing.assert_array_equal(row_sums, np.ones(4))

    # cell1/cell2 share one collapsed guide-combo column; cell3/cell4 share a
    # different one, and the two groups must NOT collide.
    _, col_idx = collapsed_x.nonzero()
    assert col_idx[0] == col_idx[1]
    assert col_idx[2] == col_idx[3]
    assert col_idx[0] != col_idx[2]


def test_collapse_guides_layer_matches_x():
    mdata = _make_dual_guide_mudata()
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        input_path = f"{d}/input.h5mu"
        output_path = f"{d}/output.h5mu"
        mdata.write(input_path)
        result = cg.collapse_guides(input_path, output_path)

    x = result["guide"].X
    layer = result["guide"].layers["guide_assignment"]
    assert sparse.issparse(layer)
    assert layer.dtype == np.uint16
    np.testing.assert_array_equal(x.toarray(), layer.toarray())


def test_collapse_guides_survives_write_read_roundtrip(tmp_path):
    mdata = _make_dual_guide_mudata()
    input_path = tmp_path / "input.h5mu"
    output_path = tmp_path / "output.h5mu"
    mdata.write(input_path)

    cg.collapse_guides(str(input_path), str(output_path))
    reloaded = mudata.read_h5mu(output_path)

    reloaded_x = reloaded["guide"].X
    assert sparse.issparse(reloaded_x)
    assert reloaded_x.dtype == np.uint16
    row_sums = np.asarray(reloaded_x.sum(axis=1)).ravel()
    np.testing.assert_array_equal(row_sums, np.ones(4))
