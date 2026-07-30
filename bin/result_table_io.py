#!/usr/bin/env python
"""Read and write inference result tables in TSV or Parquet format."""

from pathlib import Path

import pandas as pd


SUPPORTED_RESULT_FORMATS = {"tsv.gz", "parquet"}


def normalize_result_format(value):
    normalized = str(value).strip().lower().lstrip(".")
    aliases = {
        "tsv": "tsv.gz",
        "gz": "tsv.gz",
        "tsv.gz": "tsv.gz",
        "parquet": "parquet",
        "pq": "parquet",
    }
    try:
        return aliases[normalized]
    except KeyError as exc:
        choices = ", ".join(sorted(SUPPORTED_RESULT_FORMATS))
        raise ValueError(
            f"Unsupported result format {value!r}; expected one of: {choices}"
        ) from exc


def result_suffix(value):
    return ".parquet" if normalize_result_format(value) == "parquet" else ".tsv.gz"


def detect_result_format(path):
    path_string = str(path).lower()
    if path_string.endswith((".parquet", ".pq")):
        return "parquet"
    if path_string.endswith((".tsv.gz", ".tsv", ".txt.gz", ".txt")):
        return "tsv.gz"
    raise ValueError(
        f"Cannot infer result-table format from {path!r}; "
        "use a .parquet, .pq, .tsv, or .tsv.gz suffix"
    )


def read_result_table(path):
    result_format = detect_result_format(path)
    if result_format == "parquet":
        return pd.read_parquet(path, engine="pyarrow")
    return pd.read_csv(path, sep="\t")


def make_parquet_safe(frame):
    """Normalize scalar object columns so PyArrow sees one stable type."""
    normalized = frame.copy()
    for column in normalized.columns:
        if pd.api.types.is_object_dtype(normalized[column].dtype):
            normalized[column] = normalized[column].astype("string")
    return normalized


def write_result_table(frame, path, parquet_compression="zstd"):
    destination = Path(path)
    result_format = detect_result_format(destination)
    if result_format == "parquet":
        make_parquet_safe(frame).to_parquet(
            destination,
            engine="pyarrow",
            compression=parquet_compression,
            index=False,
        )
    else:
        compression = "gzip" if str(destination).lower().endswith(".gz") else None
        frame.to_csv(
            destination,
            index=False,
            sep="\t",
            compression=compression,
        )
    return destination
