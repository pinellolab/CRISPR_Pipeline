#!/usr/bin/env python3
"""Read only the parts of an inference MuData that the QC scripts use.

The inference MuData keeps the full result tables in ``uns``, and they dominate
the file: on the TAP-seq chr8 screen ``uns`` is 6.35 GB of a 6.41 GB file --
4.99 GB of that the 52,006,760-row guide table -- against 0.06 GB for the
modalities themselves. ``mudata.read_h5mu`` loads ``uns`` eagerly (``backed``
defers only ``.X``), so a script that wants ``.var`` pays for all of it.

``additional_qc_plots`` runs five such scripts over the same file, which is how
that one task read 43.3 GB and peaked at 41.9 GB of RSS to write 0.6 GB of
plots. Three of the five never touch a result table, and the two that do need a
handful of columns out of twenty-odd.

Nothing here changes what the QC scripts compute; it changes what they read.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np
import pandas as pd

# The result tables the QC scripts look for, in the order they prefer them.
RESULTS_KEY_CANDIDATES = (
    "global_analysis_per_guide_results",
    "local_analysis_per_guide_results",
    "trans_per_guide_results",
    "per_guide_results",
    "cis_per_guide_results",
    "trans_test_results",
    "test_results",
)

# The columns intended_target.py and trans.py read, plus every metric name their
# own auto-detection considers, so resolution still sees what it would have seen.
# Columns those scripts compute (p_value_adj, significant, negLog10p) are not
# read and so are not listed. Absent names are skipped, not an error.
QC_RESULT_COLUMNS = (
    "gene_id",
    "guide_id",
    "gene_name",
    "intended_target_name",
    "target_name",
    "targeting",
    "is_validated_pair",
    "log2_fc",
    "perturbo_log2_fc",
    "sceptre_log2_fc",
    "p_value",
    "perturbo_p_value",
    "sceptre_p_value",
    "q_value",
)


def _decode(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in "SO":
        return np.array(
            [v.decode() if isinstance(v, bytes) else v for v in values], dtype=object
        )
    return values


def result_keys(path: str | Path) -> list[str]:
    """Names of the ``uns`` result tables, reading none of them."""
    with h5py.File(path, "r") as handle:
        if "uns" not in handle:
            return []
        return sorted(handle["uns"].keys())


def result_table_columns(path: str | Path, key: str) -> list[str]:
    """Column names of one ``uns`` result table, reading no column."""
    with h5py.File(path, "r") as handle:
        group = handle.get(f"uns/{key}")
        if group is None:
            return []
        order = group.attrs.get("column-order")
        if order is not None:
            return [c.decode() if isinstance(c, bytes) else str(c) for c in order]
        return [name for name in group.keys() if name != "_index"]


def read_result_columns(
    path: str | Path, key: str, columns: Iterable[str]
) -> pd.DataFrame:
    """Selected columns of one ``uns`` result table.

    Decodes the encodings anndata writes -- categorical as categories plus
    codes, nullable types as values plus mask -- so the result matches what
    ``read_h5mu`` would have produced for those columns. A column the table does
    not have is skipped rather than raising, which mirrors the callers' own
    tolerance for absent optional metrics.
    """
    wanted = list(dict.fromkeys(columns))
    frame = {}
    with h5py.File(path, "r") as handle:
        group = handle.get(f"uns/{key}")
        if group is None:
            raise KeyError(f"uns/{key} not found in {path}")
        for column in wanted:
            node = group.get(column)
            if node is None:
                continue
            encoding = node.attrs.get("encoding-type", b"")
            encoding = (
                encoding.decode() if isinstance(encoding, bytes) else str(encoding)
            )
            if encoding == "categorical":
                categories = _decode(node["categories"][:])
                codes = node["codes"][:]
                values = np.full(len(codes), None, dtype=object)
                present = codes >= 0
                values[present] = categories[codes[present]]
                frame[column] = values
            elif encoding.startswith("nullable"):
                values = node["values"][:]
                mask = node["mask"][:]
                if values.dtype.kind in "iu":
                    values = values.astype("float64")
                values = np.asarray(values, dtype=object if values.dtype.kind in "SO" else values.dtype)
                if values.dtype.kind == "f":
                    values[mask] = np.nan
                else:
                    values = values.astype(object)
                    values[mask] = None
                frame[column] = values
            else:
                frame[column] = _decode(node[:])
    return pd.DataFrame(frame)


def read_mudata_without_uns(path: str | Path):
    """A MuData carrying every modality but no ``uns`` result tables.

    ``mudata.read_h5ad(path, mod=...)`` reads one modality out of an ``.h5mu``
    without touching ``uns``, so this just assembles those into a MuData and
    callers keep the ``.mod`` interface they already use. The modalities are the
    0.06 GB part of the file.
    """
    import mudata as md

    with h5py.File(path, "r") as handle:
        names = list(handle["mod"])
    return md.MuData({name: md.read_h5ad(path, mod=name) for name in names})


def resolve_result_key(
    path: str | Path, requested: str, candidates: Sequence[str]
) -> str | None:
    """The result table to use, chosen from names alone."""
    available = result_keys(path)
    if requested != "auto" and requested in available:
        return requested
    if requested != "auto":
        print(
            f"Requested results key '{requested}' not found. Falling back to "
            "auto-detection.",
            flush=True,
        )
    for candidate in candidates:
        if candidate in available:
            return candidate
    return None
