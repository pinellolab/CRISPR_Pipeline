import argparse
import pathlib
import sys
import types
from types import SimpleNamespace

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

pd = pytest.importorskip("pandas")
np = pytest.importorskip("numpy")
sparse = pytest.importorskip("scipy.sparse")
pytest.importorskip("Bio")

# crispr_ambiguous_mapping is a specialized dependency only present in the
# base-editing container image; series_to_anndata (what's under test here)
# never calls into it, so stub it out rather than skipping this test
# wholesale in environments that lack it.
if "crispr_ambiguous_mapping" not in sys.modules:
    sys.modules["crispr_ambiguous_mapping"] = types.ModuleType("crispr_ambiguous_mapping")

import base_editing_mapping as bem


def _make_guide_set(tmp_path, spacers_to_ids):
    path = tmp_path / "guide_set.tsv"
    pd.DataFrame(
        {"spacer": list(spacers_to_ids.keys()), "guide_id": list(spacers_to_ids.values())}
    ).to_csv(path, sep="\t", index=False)
    return path


def test_series_to_anndata_produces_sparse_narrow_dtype(tmp_path):
    # MultiIndex (CellBarcode, protospacer) counts, mostly zero once unstacked.
    index = pd.MultiIndex.from_tuples(
        [
            ("cellA", "AAAA"),
            ("cellA", "CCCC"),
            ("cellB", "CCCC"),
        ],
        names=["CellBarcode", "protospacer"],
    )
    counts = pd.Series([3, 0, 7], index=index)

    guide_set_fn = _make_guide_set(tmp_path, {"AAAA": "guideA", "CCCC": "guideC"})
    args = SimpleNamespace(guide_set_fn=guide_set_fn)

    adata = bem.series_to_anndata(counts, args)

    assert sparse.issparse(adata.X)
    assert adata.X.dtype == np.uint16
    assert list(adata.obs_names) == ["cellA", "cellB"]

    expected = counts.unstack(level="protospacer").fillna(0).reindex(columns=adata.var_names)
    np.testing.assert_array_equal(adata.X.toarray(), expected.to_numpy())

    # guide_id mapping applied correctly, not just column order.
    mapping = dict(zip(adata.var_names, adata.var["guide_id"]))
    assert mapping["AAAA"] == "guideA"
    assert mapping["CCCC"] == "guideC"


def test_series_to_anndata_raises_on_missing_spacer_metadata(tmp_path):
    index = pd.MultiIndex.from_tuples(
        [("cellA", "GGGG")], names=["CellBarcode", "protospacer"]
    )
    counts = pd.Series([1], index=index)
    guide_set_fn = _make_guide_set(tmp_path, {"AAAA": "guideA"})
    args = SimpleNamespace(guide_set_fn=guide_set_fn)

    with pytest.raises(ValueError, match="absent from guide metadata"):
        bem.series_to_anndata(counts, args)
