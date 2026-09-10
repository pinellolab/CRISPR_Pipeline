import pathlib
import sys

import pandas as pd


BIN_DIR = pathlib.Path(__file__).resolve().parents[1] / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from inference_target_matching import direct_target_mask


def test_direct_targets_match_gene_symbols_and_ensembl_ids():
    results = pd.DataFrame(
        {
            "intended_target_name": ["FAM83A", "ENSG000001234.7", "chr8:1-20", None],
            "gene_id": ["ENSG00000147689", "ENSG000001234", "ENSG000009999", "FAM83A"],
            "gene_name": ["FAM83A", "OTHER", "INTERVAL_GENE", "FAM83A"],
        }
    )

    assert direct_target_mask(results).tolist() == [True, True, False, False]


def test_direct_target_matching_is_case_and_whitespace_tolerant():
    results = pd.DataFrame(
        {
            "intended_target_name": [" ripk2 "],
            "gene_id": ["ENSG00000104312"],
            "gene_name": ["RIPK2"],
        }
    )

    assert direct_target_mask(results).tolist() == [True]


def test_direct_target_matching_keeps_categorical_columns_compact():
    results = pd.DataFrame(
        {
            "intended_target_name": pd.Categorical(["FAM83A", "ENSG00000104312.3", None]),
            "gene_id": pd.Categorical(["ENSG00000147689", "ENSG00000104312", "ENSG00000147689"]),
            "gene_name": pd.Categorical(["FAM83A", "RIPK2", "FAM83A"]),
        }
    )

    assert direct_target_mask(results).tolist() == [True, True, False]
