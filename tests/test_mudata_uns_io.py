import pathlib
import sys

import anndata as ad
import h5py
import mudata as mu
import numpy as np
import pandas as pd
import pytest
from scipy import sparse


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from mudata_uns_io import write_uns_patch


def _make_mudata(path, n_obs=200, n_var=50):
    X = sparse.random(n_obs, n_var, density=0.1, format="csr", dtype=np.float32, random_state=0)
    gene = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"cell{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=[f"gene{i}" for i in range(n_var)]),
    )
    mdata = mu.MuData({"gene": gene})
    mdata.uns["keep_me"] = pd.DataFrame({"x": [1, 2]})
    mdata.uns["drop_me"] = pd.DataFrame({"y": [3, 4]})
    mdata.write(path, compression="gzip")
    return mdata


def _x_bytes(path):
    with h5py.File(path, "r") as f:
        return f["mod/gene/X/data"][:].copy()


def test_write_uns_patch_adds_and_deletes_without_touching_x(tmp_path):
    input_path = tmp_path / "input.h5mu"
    output_path = tmp_path / "output.h5mu"
    _make_mudata(input_path)

    x_before = _x_bytes(input_path)

    write_uns_patch(
        input_path,
        output_path,
        updates={"new_key": pd.DataFrame({"z": [5, 6]})},
        deletes=["drop_me"],
    )

    x_after = _x_bytes(output_path)
    assert np.array_equal(x_before, x_after), "X data should be byte-identical after a /uns-only patch"

    result = mu.read_h5mu(output_path)
    assert set(result.uns.keys()) == {"keep_me", "new_key"}
    pd.testing.assert_frame_equal(pd.DataFrame(result.uns["keep_me"]), pd.DataFrame({"x": [1, 2]}))
    pd.testing.assert_frame_equal(pd.DataFrame(result.uns["new_key"]), pd.DataFrame({"z": [5, 6]}))

    # Original input file must be untouched.
    original = mu.read_h5mu(input_path)
    assert set(original.uns.keys()) == {"keep_me", "drop_me"}


def test_write_uns_patch_same_path_updates_in_place(tmp_path):
    path = tmp_path / "inplace.h5mu"
    _make_mudata(path)

    write_uns_patch(path, path, updates={"new_key": pd.DataFrame({"z": [7]})})

    result = mu.read_h5mu(path)
    assert "new_key" in result.uns
    assert "keep_me" in result.uns
