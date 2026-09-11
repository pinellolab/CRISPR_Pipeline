#!/usr/bin/env python
"""Bridge the Polars LazyFrame schema API across versions.

Polars 1.0 made ``LazyFrame.schema`` eager-only and added ``collect_schema()``;
the pipeline's base container ships 0.20, where only the property exists. The
streaming merge used ``collect_schema()`` unconditionally, so every global-format
run died in mergeMudata with ``'LazyFrame' object has no attribute
'collect_schema'``.
"""

from __future__ import annotations


class _Schema(dict):
    """A mapping with the small part of the 1.0 Schema surface this code uses."""

    def names(self):
        return list(self.keys())

    def dtypes(self):
        return list(self.values())


def lazy_schema(scan):
    """The schema of a LazyFrame, on either Polars generation."""
    if hasattr(scan, "collect_schema"):
        return scan.collect_schema()
    return _Schema(scan.schema)
