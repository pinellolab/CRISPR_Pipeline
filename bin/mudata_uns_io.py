#!/usr/bin/env python
"""Patch a MuData file's top-level .uns without rewriting its assay matrices."""

from __future__ import annotations

from pathlib import Path
import shutil
from typing import Iterable

import h5py

try:
    from anndata.io import read_elem, write_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem, write_elem


def write_uns_patch(
    input_path: str | Path,
    output_path: str | Path,
    updates: dict | None = None,
    deletes: Iterable[str] = (),
) -> None:
    """Copy ``input_path`` to ``output_path``, then add/overwrite ``updates``
    and remove ``deletes`` from the file's top-level /uns group.

    Never touches /mod/*/X or /mod/*/layers, so this is far cheaper than a
    full MuData read + write when the only change is to .uns -- the copy is
    a raw byte copy (no re-parsing/re-chunking/re-compression), and the /uns
    patch only re-serializes whatever was already in .uns plus `updates`.
    """
    input_path = Path(input_path)
    output_path = Path(output_path)
    if input_path.resolve() != output_path.resolve():
        shutil.copy(input_path, output_path)

    with h5py.File(output_path, "r+") as f:
        existing = read_elem(f["uns"]) if "uns" in f else {}
        existing.update(updates or {})
        for key in deletes:
            existing.pop(key, None)
        write_elem(f, "uns", existing)
