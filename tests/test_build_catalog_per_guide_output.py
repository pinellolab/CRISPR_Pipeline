import pathlib
import sys

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
from scipy import sparse

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

import build_catalog_per_guide_output as catalog_builder


def _make_test_mudata():
    obs = pd.DataFrame(index=["cell1", "cell2", "cell3"])
    gene = ad.AnnData(
        X=np.ones((3, 2)),
        obs=obs.copy(),
        var=pd.DataFrame({"symbol": ["SYM1", "SYM2"]}, index=["GENE1", "GENE2"]),
    )
    guide_var = pd.DataFrame(
        {
            "guide_id": ["g2", "g1"],
            "spacer": ["CCCC", "AAAA"],
            "type": ["targeting", "targeting"],
            "targeting": [True, True],
            "guide_chr": ["chr2", "chr1"],
            "guide_start": [200, 100],
            "guide_end": [220, 120],
            "strand": ["-", "+"],
            "pam": ["NGG", "NGG"],
            "intended_target_name": ["elem2", "elem1"],
            "intended_target_chr": ["chr2", "chr1"],
            "intended_target_start": [190, 90],
            "intended_target_end": [230, 130],
        },
        index=["g2", "g1"],
    )
    guide_var.index.name = "guide_id"
    guide = ad.AnnData(X=np.ones((3, 2)), obs=obs.copy(), var=guide_var)
    guide.layers["guide_assignment"] = sparse.csr_matrix(
        np.array([[1, 0], [0, 1], [0, 1]], dtype=float)
    )
    return mu.MuData({"gene": gene, "guide": guide})


def _make_input_tables():
    local = pd.DataFrame(
        {
            "gene_id": ["GENE1"],
            "guide_id": ["g1"],
            "sceptre_log2_fc": [1.2],
            "sceptre_p_value": [0.01],
            "sceptre_q_value": [0.02],
            "sceptre_fc_se": [0.3],
        }
    )
    global_results = pd.DataFrame(
        {
            "gene_id": ["GENE1", "GENE2"],
            "guide_id": ["g1", "g2"],
            "log2_fc": [0.5, -1.0],
            "p_value": [0.2, 0.0],
        }
    )
    return local, global_results


def test_create_catalog_per_guide_merge_metadata_and_counts():
    catalog = catalog_builder.create_catalog_per_guide(
        *_make_input_tables(), _make_test_mudata()
    )

    assert list(catalog.columns) == catalog_builder.OUTPUT_COLUMNS
    assert len(catalog) == 2

    g1 = catalog[(catalog["guide_id"] == "g1") & (catalog["gene_id"] == "GENE1")].iloc[0]
    assert g1["guide_sequence"] == "AAAA"
    assert g1["guide_type"] == "targeting"
    assert g1["guide_strand"] == "+"
    assert g1["gene_name"] == "SYM1"
    assert g1["nPerturbedCells"] == 2
    assert np.isclose(g1["sceptre_negLog10p"], 2.0)
    assert np.isclose(g1["perturbo_negLog10p"], -np.log10(0.2))

    g2 = catalog[(catalog["guide_id"] == "g2") & (catalog["gene_id"] == "GENE2")].iloc[0]
    assert pd.isna(g2["sceptre_log2_fc"])
    assert g2["perturbo_negLog10p"] == 300.0
    assert g2["gene_name"] == "SYM2"
    assert g2["nPerturbedCells"] == 1


def test_build_catalog_per_guide_writes_gzip_tsv(tmp_path):
    local, global_results = _make_input_tables()
    mdata = _make_test_mudata()
    local_path = tmp_path / "local_analysis_per_guide_output.tsv.gz"
    global_path = tmp_path / "global_analysis_per_guide_output.tsv.gz"
    mudata_path = tmp_path / "inference_mudata.h5mu"
    output_path = tmp_path / "catalog_per_guide_output.tsv.gz"
    local.to_csv(local_path, sep="\t", index=False, compression="gzip")
    global_results.to_csv(global_path, sep="\t", index=False, compression="gzip")
    mdata.write(mudata_path)

    catalog_builder.build_catalog_per_guide_output(
        str(local_path), str(global_path), str(mudata_path), str(output_path)
    )

    assert output_path.exists()
    observed = pd.read_csv(output_path, sep="\t")
    assert list(observed.columns) == catalog_builder.OUTPUT_COLUMNS
    assert len(observed) == 2
