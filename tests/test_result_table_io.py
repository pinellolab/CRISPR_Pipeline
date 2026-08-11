import pathlib
import sys

import pandas as pd
import pytest


REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
BIN_DIR = REPO_ROOT / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from result_table_io import categoricalize_text_columns, read_result_table, write_result_table


def test_parquet_result_tables_store_all_text_columns_as_categoricals(tmp_path):
    pytest.importorskip("pyarrow")
    frame = pd.DataFrame(
        {
            "guide_id": ["g1", "g1", "g2"],
            "gene_id": pd.Series(["GENE1", "GENE2", "GENE1"], dtype="string"),
            "p_value": [0.01, 0.02, 0.03],
        }
    )

    prepared = categoricalize_text_columns(frame)
    assert str(prepared["guide_id"].dtype) == "category"
    assert str(prepared["gene_id"].dtype) == "category"

    output = tmp_path / "results.parquet"
    write_result_table(frame, output)
    observed = read_result_table(output)
    assert observed["guide_id"].dtype.name == "category"
    assert observed["gene_id"].dtype.name == "category"
