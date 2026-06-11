import pathlib
import sys

import pytest

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from barcode_keys import qualify_barcodes


def test_qualify_barcodes_uses_measurement_set_id():
    assert qualify_barcodes(
        ["AAAC", "TTTG"],
        "IGVFDS6244NAXC",
    ) == [
        "AAAC_IGVFDS6244NAXC",
        "TTTG_IGVFDS6244NAXC",
    ]


def test_qualify_barcodes_separates_same_raw_barcode_across_batches():
    first = qualify_barcodes(["AAAC"], "IGVFDS6244NAXC")
    second = qualify_barcodes(["AAAC"], "IGVFDS8721BKRO")

    assert first[0] != second[0]


def _load_anndata_modules():
    ad = pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    pd = pytest.importorskip("pandas")
    import anndata_concat

    return ad, np, pd, anndata_concat


def _adata(barcodes):
    ad, np, pd, _ = _load_anndata_modules()
    return ad.AnnData(
        X=np.ones((len(barcodes), 1)),
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=["feature"]),
    )


def test_apply_batch_suffix_records_corrected_barcode_when_requested():
    _, _, _, anndata_concat = _load_anndata_modules()
    adata = _adata(["AAAC", "TTTG"])

    anndata_concat.apply_batch_suffix(
        adata,
        batch_num="IGVFDS6244NAXC",
        barcode_suffix="IGVFDS6244NAXC",
        seen_barcodes={},
        record_corrected_barcode=True,
    )

    assert adata.obs_names.tolist() == [
        "AAAC_IGVFDS6244NAXC",
        "TTTG_IGVFDS6244NAXC",
    ]
    assert adata.obs["corrected_barcode"].tolist() == ["AAAC", "TTTG"]


def test_apply_batch_suffix_rejects_duplicate_batch_qualified_barcode():
    _, _, _, anndata_concat = _load_anndata_modules()
    seen_barcodes = {}
    first = _adata(["AAAC"])
    second = _adata(["AAAC"])

    anndata_concat.apply_batch_suffix(
        first,
        batch_num="batch_1",
        barcode_suffix="batch_1",
        seen_barcodes=seen_barcodes,
    )

    with pytest.raises(ValueError, match="duplicated across mapping outputs"):
        anndata_concat.apply_batch_suffix(
            second,
            batch_num="batch_1",
            barcode_suffix="batch_1",
            seen_barcodes=seen_barcodes,
        )
