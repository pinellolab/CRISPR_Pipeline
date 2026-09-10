#!/usr/bin/env python3
"""Identifier-aware matching for guide intended-target inference results."""

from __future__ import annotations

import re

import numpy as np
import pandas as pd


_ENSEMBL_VERSION = re.compile(r"^(ENS[A-Z]*G\d+)\.\d+$", re.IGNORECASE)


def normalize_target_identifier(values: pd.Series) -> pd.Series:
    """Normalize gene symbols/IDs without treating missing values as strings."""
    normalized = values.astype("string").str.strip().str.upper()
    return normalized.str.replace(_ENSEMBL_VERSION, r"\1", regex=True)


def _normalize_scalar(value: object) -> str | None:
    if pd.isna(value):
        return None
    normalized = str(value).strip().upper()
    match = _ENSEMBL_VERSION.match(normalized)
    return match.group(1) if match else normalized


def _identifier_equals(left: pd.Series, right: pd.Series) -> pd.Series:
    """Compare identifiers without expanding large categorical columns to strings."""
    if isinstance(left.dtype, pd.CategoricalDtype) and isinstance(
        right.dtype, pd.CategoricalDtype
    ):
        right_codes_by_value = {
            _normalize_scalar(value): code
            for code, value in enumerate(right.cat.categories)
        }
        left_to_right = np.asarray(
            [right_codes_by_value.get(_normalize_scalar(value), -2) for value in left.cat.categories],
            dtype=np.int32,
        )
        left_codes = left.cat.codes.to_numpy()
        right_codes = right.cat.codes.to_numpy()
        expected = np.full(left_codes.shape, -2, dtype=np.int32)
        valid = left_codes >= 0
        expected[valid] = left_to_right[left_codes[valid]]
        return pd.Series(valid & (expected == right_codes), index=left.index)

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
