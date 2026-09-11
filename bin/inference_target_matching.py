#!/usr/bin/env python3
"""Identifier-aware matching for guide intended-target inference results."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd


_ENSEMBL_VERSION = re.compile(r"^(ENS[A-Z]*G\d+)\.\d+$", re.IGNORECASE)


def normalize_target_identifier(values: pd.Series) -> pd.Series:
    """Normalize gene symbols/IDs without treating missing values as strings.

    An identifier that is empty once stripped carries no information, so it is
    normalized to missing rather than to the empty string.  Several upstream
    steps fill an absent name with ``""`` (``analysis_output_formatting``'s
    guide metadata and gene-name fill, and the Polars fast path's
    ``fill_null("")``), and two absent names must never be read as the same
    gene.
    """
    normalized = values.astype("string").str.strip().str.upper()
    normalized = normalized.str.replace(_ENSEMBL_VERSION, r"\1", regex=True)
    return normalized.mask(normalized.eq(""), pd.NA)


def _normalize_scalar(value: object) -> str | None:
    if pd.isna(value):
        return None
    normalized = str(value).strip().upper()
    if not normalized:
        return None
    match = _ENSEMBL_VERSION.match(normalized)
    return match.group(1) if match else normalized


def _categorical_identifier_codes(
    categories: pd.Index, vocabulary: dict[str, int], missing: int
) -> np.ndarray:
    """Map each category to a shared id, or ``missing`` when it carries none."""
    codes = np.full(len(categories), missing, dtype=np.int64)
    for position, value in enumerate(categories):
        key = _normalize_scalar(value)
        if key is None:
            continue
        codes[position] = vocabulary.setdefault(key, len(vocabulary))
    return codes


def _identifier_equals(left: pd.Series, right: pd.Series) -> pd.Series:
    """Compare identifiers without expanding large categorical columns to strings."""
    if isinstance(left.dtype, pd.CategoricalDtype) and isinstance(
        right.dtype, pd.CategoricalDtype
    ):
        # Both sides are resolved against one shared vocabulary of normalized
        # identifiers, so two categories that normalize to the same gene (say
        # ``ENSG1`` and ``ENSG1.3``) compare equal here exactly as they do on
        # the string path below.
        vocabulary: dict[str, int] = {}
        left_ids = _categorical_identifier_codes(left.cat.categories, vocabulary, -1)
        right_ids = _categorical_identifier_codes(right.cat.categories, vocabulary, -2)
        left_codes = left.cat.codes.to_numpy()
        right_codes = right.cat.codes.to_numpy()
        # A trailing sentinel lets an unset code (-1) index the lookup safely,
        # including when a column has no categories at all.
        left_lookup = np.append(left_ids, -1)
        right_lookup = np.append(right_ids, -2)
        left_values = left_lookup[np.where(left_codes >= 0, left_codes, len(left_ids))]
        right_values = right_lookup[
            np.where(right_codes >= 0, right_codes, len(right_ids))
        ]
        return pd.Series(
            (left_values >= 0) & (left_values == right_values), index=left.index
        )

    left_normalized = normalize_target_identifier(left)
    right_normalized = normalize_target_identifier(right)
    return left_normalized.notna() & right_normalized.notna() & left_normalized.eq(right_normalized)


def direct_target_mask(
    results: pd.DataFrame,
    target_col: str = "intended_target_name",
    gene_id_col: str = "gene_id",
    gene_name_col: str = "gene_name",
) -> pd.Series:
    """Return rows where an intended target matches gene ID or gene symbol.

    Inference output commonly stores Ensembl IDs in ``gene_id`` while guide
    designs store symbols (for example, ``FAM83A``) in
    ``intended_target_name``.  Genomic-element names remain unmatched unless
    they genuinely equal a tested gene identifier.
    """
    missing = [column for column in (target_col, gene_id_col) if column not in results]
    if missing:
        raise KeyError(f"Missing required direct-target columns: {', '.join(missing)}")

    target = results[target_col]
    matched = _identifier_equals(target, results[gene_id_col])

    if gene_name_col in results:
        matched |= _identifier_equals(target, results[gene_name_col])

    return matched.fillna(False)
