import pathlib
import sys

import pandas as pd


BIN_DIR = pathlib.Path(__file__).resolve().parents[1] / "bin"
if str(BIN_DIR) not in sys.path:
    sys.path.insert(0, str(BIN_DIR))

from streaming_catalog_io import try_write_enriched_parquet_catalog


def test_streaming_catalog_normalizes_categorical_and_string_join_keys(tmp_path):
    local = pd.DataFrame(
        {
            "gene_id": pd.Categorical(["gene1", "gene2"]),
            "guide_id": pd.Categorical(["guide1", "guide2"]),
            "local_score": [1.0, 2.0],
        }
    )
    global_results = pd.DataFrame(
        {
            "gene_id": pd.Series(["gene1", "gene2"], dtype="string"),
            "guide_id": pd.Series(["guide1", "guide2"], dtype="string"),
            "global_score": [3.0, 4.0],
        }
    )
    local_path = tmp_path / "local.parquet"
    global_path = tmp_path / "global.parquet"
    output_path = tmp_path / "catalog.parquet"
    local.to_parquet(local_path, index=False)
    global_results.to_parquet(global_path, index=False)

    used_fast_path = try_write_enriched_parquet_catalog(
        local_path=local_path,
        global_path=global_path,
        output_path=output_path,
        join_columns=["gene_id", "guide_id"],
        local_metric_columns=["local_score"],
        global_required_columns=["global_score"],
        output_columns=[
            "gene_id",
            "guide_id",
            "local_score",
            "global_score",
        ],
        sort_columns=["gene_id", "guide_id"],
    )

    assert used_fast_path
    observed = pd.read_parquet(output_path)
    assert observed["local_score"].tolist() == [1.0, 2.0]
    assert observed["global_score"].tolist() == [3.0, 4.0]
