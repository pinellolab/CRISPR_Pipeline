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


def test_normalize_feature_index_promotes_legacy_flash_guide_ids():
    ad, np, pd, anndata_concat = _load_anndata_modules()
    adata = ad.AnnData(
        X=np.ones((1, 2)),
        obs=pd.DataFrame(index=["AAAC"]),
        var=pd.DataFrame(
            {"guide_id": ["guide-1", "guide-2"]},
            index=["AACCGGTTAACCGGTTAACC", "TTGGCCAATTGGCCAATTGG"],
        ),
    )

    anndata_concat.normalize_feature_index(adata, "legacy FLASH output")

    assert adata.var_names.tolist() == ["guide-1", "guide-2"]
    assert adata.var_names.name == "guide_id"
    assert "guide_id" not in adata.var.columns
    assert adata.var["spacer"].tolist() == [
        "AACCGGTTAACCGGTTAACC",
        "TTGGCCAATTGGCCAATTGG",
    ]


def test_flash_mapper_writes_guide_ids_as_feature_index():
    ad, np, pd, _ = _load_anndata_modules()
    import flash_base_editing_mapping

    guide_df = pd.DataFrame({"guide_id": ["guide-1", "guide-2"]})
    adata = flash_base_editing_mapping.build_output_anndata(
        counts_by_cell={1: {0: {11, 12}, 1: {13}}},
        barcode_key_to_str={1: "AAAC"},
        guide_df=guide_df,
        raw_spacers=["AACCGGTTAACCGGTTAACC", "TTGGCCAATTGGCCAATTGG"],
        output_to_unique=np.asarray([0, 1], dtype=np.int32),
    )

    assert isinstance(adata, ad.AnnData)
    assert adata.var_names.tolist() == ["guide-1", "guide-2"]
    assert adata.var_names.name == "guide_id"
    assert adata.var["spacer"].tolist() == [
        "AACCGGTTAACCGGTTAACC",
        "TTGGCCAATTGGCCAATTGG",
    ]
    assert adata.X.toarray().tolist() == [[2, 1]]


def test_legacy_flash_concat_keeps_guide_id_for_metadata_merge(tmp_path):
    ad, np, pd, anndata_concat = _load_anndata_modules()
    sparse = pytest.importorskip("scipy.sparse")
    input_paths = []
    for batch in ("A", "B"):
        adata = ad.AnnData(
            X=sparse.csr_matrix(np.ones((1, 2))),
            obs=pd.DataFrame(index=[f"cell-{batch}"]),
            var=pd.DataFrame(
                {"guide_id": ["guide-1", "guide-2"]},
                index=["AACCGGTTAACCGGTTAACC", "TTGGCCAATTGGCCAATTGG"],
            ),
        )
        anndata_concat.normalize_feature_index(adata, f"legacy FLASH batch {batch}")
        input_path = tmp_path / f"{batch}.h5ad"
        adata.write_h5ad(input_path)
        input_paths.append(input_path)

    output_path = tmp_path / "concatenated.h5ad"
    ad.experimental.concat_on_disk(
        input_paths,
        join="outer",
        index_unique=None,
        out_file=output_path,
    )
    concatenated = ad.read_h5ad(output_path)
    # anndata_concat restores the name captured from the first normalized input.
    concatenated.var_names.name = "guide_id"

    guide_var = concatenated.var.reset_index()
    metadata = pd.DataFrame(
        {
            "guide_id": ["guide-1", "guide-2"],
            "targeting": [True, False],
        }
    )
    merged = guide_var.merge(metadata, on="guide_id", validate="one_to_one")

    assert concatenated.var_names.name == "guide_id"
    assert concatenated.var_names.tolist() == ["guide-1", "guide-2"]
    assert merged["targeting"].tolist() == [True, False]
