"""Bounded conversion keeps one stable public Arrow schema across batches."""
from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import mudata as mu
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
import perturbo_v2_pipeline_adapter as adapter


def _prepared(path: Path) -> pd.DataFrame:
    obs = pd.DataFrame(index=["cell"])
    gene = ad.AnnData(X=np.zeros((1, 1)), obs=obs.copy(), var=pd.DataFrame(index=["G1"]))
    guide_var = pd.DataFrame(
        {
            "intended_target_key": ["control", "target"],
            "intended_target_name": ["control", "Target A"],
            "intended_target_chr": [None, "chr1"],
            "intended_target_start": pd.Series([pd.NA, 101], dtype="Int64"),
            "intended_target_end": pd.Series([pd.NA, 120], dtype="Int64"),
        },
        index=["nt", "g1"],
    )
    guide = ad.AnnData(X=np.zeros((1, 2)), obs=obs.copy(), var=guide_var)
    mu.MuData({"gene": gene, "guide": guide}).write_h5mu(path)
    return guide_var.reset_index(drop=True)


def _raw_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "element": ["control", "target"],
            "gene": ["G1", "G2"],
            "posterior_mean": np.array([0.0, np.log(2.0)], dtype=np.float32),
            "posterior_scale": np.array([0.25, 0.5], dtype=np.float32),
            "posterior_prob": np.array([0.9, 0.8], dtype=np.float64),
            "crt_saddlepoint_p_value": np.array([np.nan, 0.01], dtype=np.float64),
        }
    )


def _write(raw: Path, output: Path, prepared: Path, lookup: pd.DataFrame) -> None:
    scratch = output.parent / "scratch"
    scratch.mkdir(exist_ok=True)
    adapter._write_bounded_effects(
        raw,
        output,
        inference_type="element",
        prepared_path=None,
        target_lookup=adapter.get_target_lookup(lookup),
        guide_name_map={},
        crt=True,
        scratch_dir=scratch,
        max_bh_working_bytes=1 << 20,
        compact_floats=False,
        batch_rows=1,
    )


def test_schema_survives_all_null_first_metadata_batch_and_category_changes(tmp_path):
    prepared = tmp_path / "prepared.h5mu"
    lookup = _prepared(prepared)
    raw_frame = _raw_frame()
    raw = tmp_path / "raw.parquet"
    raw_table = pa.Table.from_pandas(raw_frame, preserve_index=False)
    with pq.ParquetWriter(raw, raw_table.schema) as writer:
        writer.write_table(raw_table.slice(0, 0))
        writer.write_table(raw_table.slice(0, 1))
        writer.write_table(raw_table.slice(1, 1))
    assert pq.ParquetFile(raw).metadata.row_group(0).num_rows == 0
    output = tmp_path / "converted.parquet"

    _write(raw, output, prepared, lookup)

    schema = pq.ParquetFile(output).schema_arrow
    for name in ("gene_id", "intended_target_name", "intended_target_chr"):
        assert pa.types.is_dictionary(schema.field(name).type)
        assert pa.types.is_string(schema.field(name).type.value_type)
    assert schema.field("intended_target_start").type == pa.int64()
    assert schema.field("intended_target_end").type == pa.int64()
    assert schema.field("log2_fc").type == pa.float32()
    assert schema.field("perturbo_fc_se").type == pa.float32()
    assert schema.field("p_value").type == pa.float64()
    assert schema.field("perturbo_q_value").type == pa.float64()

    observed = pd.read_parquet(output)
    expected = adapter._convert_common_effect_columns(raw_frame, crt=True).rename(
        columns={"element": "intended_target_key"}
    ).merge(adapter.get_target_lookup(lookup), on="intended_target_key", how="left")
    expected["perturbo_q_value"] = adapter._bh_adjust(expected["p_value"])
    expected = expected[
        [
            "gene_id", "intended_target_name", "intended_target_chr",
            "intended_target_start", "intended_target_end", "log2_fc",
            "perturbo_fc_se", "p_value", "perturbo_posterior_prob",
            "perturbo_q_value",
        ]
    ]
    for column in ("gene_id", "intended_target_name", "intended_target_chr"):
        pd.testing.assert_series_equal(
            observed[column].astype("string"), expected[column].astype("string"),
            check_names=False,
        )
    for column in expected.columns.difference(
        ["gene_id", "intended_target_name", "intended_target_chr"]
    ):
        np.testing.assert_allclose(
            pd.to_numeric(observed[column], errors="coerce"),
            pd.to_numeric(expected[column], errors="coerce"),
            rtol=1e-6,
            equal_nan=True,
        )
    assert pd.isna(observed.loc[0, "intended_target_chr"])
    assert observed.loc[1, "intended_target_chr"] == "chr1"


def test_empty_family_keeps_the_same_declared_schema(tmp_path):
    prepared = tmp_path / "prepared.h5mu"
    lookup = _prepared(prepared)
    raw = tmp_path / "empty.parquet"
    _raw_frame().iloc[:0].to_parquet(raw, index=False)
    output = tmp_path / "converted.parquet"

    _write(raw, output, prepared, lookup)

    schema = pq.ParquetFile(output).schema_arrow
    assert pq.ParquetFile(output).metadata.num_rows == 0
    assert pa.types.is_dictionary(schema.field("gene_id").type)
    assert pa.types.is_dictionary(schema.field("intended_target_chr").type)
    assert schema.field("intended_target_start").type == pa.int64()
    assert schema.field("log2_fc").type == pa.float32()
    assert schema.field("p_value").type == pa.float64()


def test_recovery_rejects_an_unmanifested_requested_table(tmp_path):
    import json
    from bounded_perturbo_results import parquet_manifest
    import perturbo_v2_pipeline_adapter as adapter

    raw = tmp_path / "element_effects.parquet"
    pd.DataFrame({"p_value": [0.1]}).to_parquet(raw, index=False)
    (tmp_path / "conversion_manifest.json").write_text(json.dumps({
        "raw_parquets": [parquet_manifest(raw)]
    }))
    pd.DataFrame({"p_value": [0.2]}).to_parquet(
        tmp_path / "element_effects_requested_pairs.parquet", index=False
    )
    import pytest
    with pytest.raises(ValueError, match="Unverified requested-pairs"):
        adapter._verify_durable_fit(tmp_path)
