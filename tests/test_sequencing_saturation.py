import importlib.util
import struct
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


SCRIPT = Path(__file__).parents[1] / "bin" / "sequencing_saturation.py"
SPEC = importlib.util.spec_from_file_location("sequencing_saturation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_bus(path, records, barcode_length=4, umi_length=4):
    text = b"synthetic BUS"
    with open(path, "wb") as handle:
        handle.write(b"BUS\x00")
        handle.write(struct.pack("<IIII", 1, barcode_length, umi_length, len(text)))
        handle.write(text)
        for record in records:
            handle.write(struct.pack("<QQiIII", *record, 0, 0))


def test_tenx_style_saturation_endpoint(tmp_path):
    mapping = tmp_path / "B1_ks_transcripts_out"
    mapping.mkdir()
    (mapping / "transcripts.txt").write_text("tx0\ntx1\n")
    (mapping / "matrix.ec").write_text("0\t0\n1\t1\n")
    t2g = tmp_path / "t2g.txt"
    t2g.write_text("tx0\tgene0\ntx1\tgene1\n")
    b0 = MODULE.encode_dna("AAAA")
    b1 = MODULE.encode_dna("AAAC")
    _write_bus(
        mapping / "output.unfiltered.bus",
        [
            (b0, 1, 0, 3),
            (b0, 2, 1, 1),
            (b1, 1, 0, 2),
        ],
    )

    ec_map, n_genes = MODULE.ec_gene_map(mapping, t2g)
    assert n_genes == 2
    result = MODULE.analyze_bus(
        mapping / "output.unfiltered.bus",
        {b0, b1},
        ec_map,
        np.array([0.5, 1.0]),
        chunk_records=2,
    )
    assert result["cells"] == 2
    assert result["usable_reads"] == 6
    assert result["expected_umis"][-1] == 3
    assert result["saturation"][-1] == 0.5
    assert result["median_umis"][-1] == 1.5
    assert result["median_genes"][-1] == 1.5


def test_filtered_barcode_suffix_recovery(tmp_path):
    filtered = ad.AnnData(
        X=np.ones((2, 1)),
        obs=pd.DataFrame({"batch": ["B1", "B1"]}, index=["AAAA_B1", "AAAC_B1"]),
    )
    covariates = pd.DataFrame({"batch": ["B1"], "barcode_key": ["B1"]})
    assert MODULE.selected_barcodes(filtered, "B1", covariates) == {
        MODULE.encode_dna("AAAA"),
        MODULE.encode_dna("AAAC"),
    }


def test_noncontiguous_equivalence_class_ids(tmp_path):
    mapping = tmp_path / "B1_ks_transcripts_out"
    mapping.mkdir()
    (mapping / "transcripts.txt").write_text("tx0\ntx1\n")
    (mapping / "matrix.ec").write_text("0\t0\n5\t1\n")
    t2g = tmp_path / "t2g.txt"
    t2g.write_text("tx0\tgene0\ntx1\tgene1\n")

    ec_map, n_genes = MODULE.ec_gene_map(mapping, t2g)

    assert len(ec_map) == 6
    assert ec_map[0] >= 0
    assert ec_map[5] >= 0
    assert n_genes == 2
