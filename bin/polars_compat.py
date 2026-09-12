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

def _supports(method: str, name: str) -> bool:
    """Whether this Polars accepts ``name`` on ``LazyFrame.method``.

    Returns False when Polars is absent. This module is imported by scripts that
    run inside the PerTurbo image, which ships no Polars at all: resolving these
    constants must not turn an optional dependency into a required one.
    """
    import inspect

    try:
        import polars as pl
    except ImportError:
        return False

    try:
        return name in inspect.signature(getattr(pl.LazyFrame, method)).parameters
    except (TypeError, ValueError, AttributeError):
        return False


# Polars grew ``maintain_order`` on joins after 1.0; the pipeline container ships
# 0.20.31, where passing it raises TypeError. Older Polars preserves the left
# frame's order for these joins anyway, so omitting the argument there gives the
# same result rather than merely a similar one. Spread these at the call site:
#     frame.join(other, on=..., how="left", **JOIN_ORDER_LEFT)
JOIN_ORDER_LEFT = {"maintain_order": "left"} if _supports("join", "maintain_order") else {}
UNIQUE_ORDER = {"maintain_order": True} if _supports("unique", "maintain_order") else {}
