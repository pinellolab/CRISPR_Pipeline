#!/usr/bin/env python3

from __future__ import annotations

import numpy as np
import pandas as pd

CONTROL_TYPES = {"safe-targeting", "non-targeting", "negative control"}
NON_TARGETING_PREFIX = "non-targeting|"


def _to_bool_series(values: pd.Series) -> pd.Series:
    """Normalize targeting-like values to bool."""
    if values is None:
        return pd.Series(dtype=bool)
    return values.map(
        lambda x: (
            bool(x)
            if isinstance(x, (bool, np.bool_))
            else str(x).strip().lower() in {"true", "1", "t", "yes"}
        )
        if not pd.isna(x)
        else False
    )


def _normalize_numeric(value):
    if pd.isna(value):
        return pd.NA
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return pd.NA


def _string_or_na(value):
    if pd.isna(value):
        return pd.NA
    return str(value)


def _build_target_key(name, chrom, start, end) -> str:
    has_coords = not (pd.isna(chrom) or pd.isna(start) or pd.isna(end))
    if has_coords:
        return f"coords::{name}::{chrom}::{int(start)}::{int(end)}"
    return f"name_only::{name}"


def _control_mask(df: pd.DataFrame) -> pd.Series:
    targeting = _to_bool_series(df.get("targeting", pd.Series(index=df.index, data=False)))
    guide_type = (
        df.get("type", pd.Series(index=df.index, data="", dtype=object))
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )
    names = (
        df.get("intended_target_name", pd.Series(index=df.index, data="", dtype=object))
        .fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )
    return (
        (~targeting)
        | guide_type.isin(CONTROL_TYPES)
        | names.eq("non-targeting")
        | names.str.startswith(NON_TARGETING_PREFIX)
    )


def _median_guides_per_targeting_element(df: pd.DataFrame, control_mask: pd.Series) -> int:
    targeting_df = df.loc[~control_mask].copy()
    if targeting_df.empty:
        return 1

    targeting_df["_tmp_name"] = targeting_df["intended_target_name"].map(_string_or_na)
    targeting_df["_tmp_chr"] = targeting_df["intended_target_chr"].map(_string_or_na)
    targeting_df["_tmp_start"] = targeting_df["intended_target_start"].map(_normalize_numeric)
    targeting_df["_tmp_end"] = targeting_df["intended_target_end"].map(_normalize_numeric)

    targeting_df["_tmp_key"] = targeting_df.apply(
        lambda row: _build_target_key(
            row["_tmp_name"],
            row["_tmp_chr"],
            row["_tmp_start"],
            row["_tmp_end"],
        ),
        axis=1,
    )

    counts = targeting_df.groupby("_tmp_key", observed=True).size()
    if counts.empty:
        return 1

    median_count = int(np.median(counts.values))
    return max(1, median_count)


def annotate_intended_target_groups(
    guide_var: pd.DataFrame,
    non_targeting_label: str = "non-targeting",
) -> pd.DataFrame:
    """
    Canonically annotate guide metadata with coordinate-aware target grouping.

    Output columns:
      - intended_target_name
      - intended_target_chr
      - intended_target_start
      - intended_target_end
      - intended_target_key

    Rules:
      - targeting guides: key = name+chr+start+end when all coordinates exist,
        else name-only fallback key.
      - control guides: grouped as non-targeting|N where group size equals median
        guides-per-targeting-element, with deterministic ordering by guide_id.
      - no guide may remain assigned to exact 'non-targeting'.
    """
    df = guide_var.copy()

    if "guide_id" not in df.columns:
        df["guide_id"] = df.index.astype(str)

    if "intended_target_name" not in df.columns:
        targeting = _to_bool_series(df.get("targeting", pd.Series(index=df.index, data=False)))
        df["intended_target_name"] = np.where(targeting, df["guide_id"].astype(str), non_targeting_label)

    if "intended_target_chr" not in df.columns:
        df["intended_target_chr"] = pd.NA
    if "intended_target_start" not in df.columns:
        df["intended_target_start"] = pd.NA
    if "intended_target_end" not in df.columns:
        df["intended_target_end"] = pd.NA

    df["guide_id"] = df["guide_id"].astype(str)
    df["intended_target_name"] = df["intended_target_name"].map(_string_or_na)
    df["intended_target_chr"] = df["intended_target_chr"].map(_string_or_na)
    df["intended_target_start"] = df["intended_target_start"].map(_normalize_numeric).astype("Int64")
    df["intended_target_end"] = df["intended_target_end"].map(_normalize_numeric).astype("Int64")

    control_mask = _control_mask(df)
    n_controls = int(control_mask.sum())

    if n_controls > 0:
        median_size = _median_guides_per_targeting_element(df, control_mask)

        control_idx_sorted = (
            df.loc[control_mask, ["guide_id"]]
            .sort_values("guide_id", kind="stable")
            .index
        )
        group_numbers = np.arange(n_controls) // median_size + 1

        control_names = [f"{non_targeting_label}|{group}" for group in group_numbers]
        control_chrom = [non_targeting_label] * n_controls
        control_start = pd.Series(group_numbers, index=control_idx_sorted, dtype="Int64")
        control_end = pd.Series(group_numbers, index=control_idx_sorted, dtype="Int64")

        df.loc[control_idx_sorted, "intended_target_name"] = control_names
        df.loc[control_idx_sorted, "intended_target_chr"] = control_chrom
        df.loc[control_idx_sorted, "intended_target_start"] = control_start
        df.loc[control_idx_sorted, "intended_target_end"] = control_end

    if (df["intended_target_name"] == non_targeting_label).any():
        raise ValueError(
            "Invalid intended_target_name assignment: found exact 'non-targeting' after bucketing."
        )

    df["intended_target_key"] = df.apply(
        lambda row: _build_target_key(
            row["intended_target_name"],
            row["intended_target_chr"],
            row["intended_target_start"],
            row["intended_target_end"],
        ),
        axis=1,
    )

    if df["intended_target_key"].isna().any():
        raise ValueError("intended_target_key contains NA values after annotation.")

    return df


def _ensure_pairs_required_columns(pairs_df: pd.DataFrame) -> pd.DataFrame:
    out = pairs_df.copy()
    for col in [
        "intended_target_name",
        "intended_target_chr",
        "intended_target_start",
        "intended_target_end",
        "intended_target_key",
    ]:
        if col not in out.columns:
            out[col] = pd.NA
    return out


def enrich_pairs_with_target_metadata(
    pairs_df: pd.DataFrame,
    guide_var: pd.DataFrame,
) -> pd.DataFrame:
    """Fill/normalize target metadata columns in a pairs table from guide metadata."""
    if "guide_id" not in pairs_df.columns:
        raise ValueError("pairs_to_test must include guide_id column")

    lookup = guide_var.reset_index(drop=True)
    lookup = lookup[
        [
            "guide_id",
            "intended_target_name",
            "intended_target_chr",
            "intended_target_start",
            "intended_target_end",
            "intended_target_key",
        ]
    ].drop_duplicates(subset=["guide_id"], keep="first")

    out = _ensure_pairs_required_columns(pairs_df)
    out["guide_id"] = out["guide_id"].astype(str)

    merged = out.merge(
        lookup,
        on="guide_id",
        how="left",
        suffixes=("", "_from_guide"),
    )

    for col in [
        "intended_target_name",
        "intended_target_chr",
        "intended_target_start",
        "intended_target_end",
    ]:
        merged[col] = merged[col].where(merged[col].notna(), merged[f"{col}_from_guide"])

    merged["intended_target_name"] = merged["intended_target_name"].map(_string_or_na)
    merged["intended_target_chr"] = merged["intended_target_chr"].map(_string_or_na)
    merged["intended_target_start"] = merged["intended_target_start"].map(_normalize_numeric).astype("Int64")
    merged["intended_target_end"] = merged["intended_target_end"].map(_normalize_numeric).astype("Int64")

    merged["intended_target_key"] = merged.apply(
        lambda row: _build_target_key(
            row["intended_target_name"],
            row["intended_target_chr"],
            row["intended_target_start"],
            row["intended_target_end"],
        ),
        axis=1,
    )

    merged = merged.drop(
        columns=[
            "intended_target_name_from_guide",
            "intended_target_chr_from_guide",
            "intended_target_start_from_guide",
            "intended_target_end_from_guide",
            "intended_target_key_from_guide",
        ],
        errors="ignore",
    )

    if merged["intended_target_key"].isna().any():
        missing = merged.loc[merged["intended_target_key"].isna(), "guide_id"].drop_duplicates().tolist()
        raise ValueError(
            "Unable to populate intended_target_key for pairs_to_test guides: "
            + ", ".join(map(str, missing[:20]))
        )

    return merged


def get_target_lookup(guide_var: pd.DataFrame) -> pd.DataFrame:
    """Return unique intended_target_key -> display metadata mapping."""
    required = [
        "intended_target_key",
        "intended_target_name",
        "intended_target_chr",
        "intended_target_start",
        "intended_target_end",
    ]
    missing = [col for col in required if col not in guide_var.columns]
    if missing:
        raise ValueError(f"guide metadata missing required columns: {', '.join(missing)}")

    lookup = (
        guide_var[required]
        .drop_duplicates(subset=["intended_target_key"], keep="first")
        .reset_index(drop=True)
    )

    if lookup["intended_target_key"].isna().any():
        raise ValueError("guide metadata contains NA intended_target_key values")

    return lookup
