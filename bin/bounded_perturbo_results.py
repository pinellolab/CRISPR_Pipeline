#!/usr/bin/env python3
"""Bounded-memory primitives for PerTurbo result conversion.

This module deliberately keeps p- and q-values as float64.  It does not decide
which rows belong to a hypothesis family; callers must pass one parquet file
containing exactly one family.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


def sha256_file(path: str | Path, *, block_size: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def parquet_manifest(path: str | Path) -> dict:
    """Describe a durable raw parquet input sufficiently for fail-closed resume."""
    import pyarrow.parquet as pq

    source = Path(path)
    parquet = pq.ParquetFile(source)
    return {
        "path": source.name,
        "size_bytes": source.stat().st_size,
        "sha256": sha256_file(source),
        "rows": parquet.metadata.num_rows,
        "row_groups": parquet.metadata.num_row_groups,
        "schema": str(parquet.schema_arrow),
    }


def write_manifest(path: str | Path, payload: dict) -> None:
    """Atomically write provenance; a crash never leaves a valid-looking file."""
    destination = Path(path)
    temporary = destination.with_name(destination.name + ".partial")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, destination)


def verify_manifest(directory: str | Path, payload: dict) -> None:
    """Fail closed if a conversion-only resume does not see identical raw data."""
    directory = Path(directory)
    for expected in payload["raw_parquets"]:
        observed = parquet_manifest(directory / expected["path"])
        if observed != expected:
            raise ValueError(
                f"Raw PerTurbo artifact failed provenance verification: "
                f"{expected['path']}"
            )


def write_bh_sidecar(
    input_path: str | Path,
    output_path: str | Path,
    *,
    p_column: str = "p_value",
    row_id_column: str = "_family_row_id",
    q_column: str = "perturbo_q_value",
    max_working_bytes: int = 24 << 30,
) -> Path:
    """Write row-id/q pairs for one complete BH family.

    Only numeric arrays are materialized.  The conservative budget accounts for
    the p-values, valid-row indices, sort permutation, and q-values.  Refusing a
    family that exceeds the budget is preferable to silently loading wide result
    frames or computing BH independently per batch.
    """
    import numpy as np
    import pyarrow as pa
    import pyarrow.parquet as pq

    source = Path(input_path)
    destination = Path(output_path)
    temporary = destination.with_name(destination.name + ".partial")
    parquet = pq.ParquetFile(source)
    n = parquet.metadata.num_rows
    # p, valid row ids, stable sort order, ranked p/q, final q, and bounded
    # vectorized temporaries. This intentionally overstates the steady state.
    estimated = n * 64
    if estimated > max_working_bytes:
        raise MemoryError(
            f"BH family has an estimated {estimated:,} numeric working bytes, above "
            f"the configured {max_working_bytes:,}-byte budget"
        )
    p = np.empty(n, dtype=np.float64)
    offset = 0
    for batch in parquet.iter_batches(columns=[p_column]):
        values = batch.column(0).to_numpy(zero_copy_only=False).astype(np.float64, copy=False)
        p[offset : offset + len(values)] = values
        offset += len(values)
    if offset != n:
        raise ValueError("Parquet metadata and p-value row counts differ")
    finite = np.isfinite(p)
    # This primitive accepts only valid probabilities or missing entries. The
    # adapter clips its BH-only projection first to preserve legacy behavior.
    if np.any((p[finite] < 0.0) | (p[finite] > 1.0)):
        raise ValueError("p-values must be between 0 and 1")
    valid_rows = np.flatnonzero(finite)
    order = np.argsort(p[valid_rows], kind="mergesort")
    ranked_p = p[valid_rows[order]]
    del p, finite
    m = ranked_p.size
    ranked_q = ranked_p
    # Match scipy's multiply-by-(m/rank) order, including subnormal p-values,
    # without a second whole-family rank or scaled-p array.
    for start in range(0, m, 250_000):
        stop = min(start + 250_000, m)
        ranked_q[start:stop] *= float(m) / np.arange(start + 1, stop + 1, dtype=np.float64)
    np.minimum.accumulate(ranked_q[::-1], out=ranked_q[::-1])
    np.minimum(ranked_q, 1.0, out=ranked_q)
    q = np.full(n, np.nan, dtype=np.float64)
    q[valid_rows[order]] = ranked_q
    del valid_rows, order, ranked_q, ranked_p
    schema = pa.schema([(row_id_column, pa.int64()), (q_column, pa.float64())])
    with pq.ParquetWriter(temporary, schema, compression="zstd") as writer:
        for start in range(0, n, 250_000):
            stop = min(start + 250_000, n)
            writer.write_table(pa.table({
                row_id_column: np.arange(start, stop, dtype=np.int64),
                q_column: q[start:stop],
            }, schema=schema))
    os.replace(temporary, destination)
    return destination
