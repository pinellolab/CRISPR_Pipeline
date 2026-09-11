#!/usr/bin/env python
"""The cell-level covariates both inference methods condition on.

They have to be agreed in one place and written where both methods look, or the
comparison between them is not a comparison of methods. PerTurbo reads named
columns from the analysed modality's ``obs``; SCEPTRE reads the MuData's top-level
``obs``, which its R reader exposes as ``colData`` and feeds to
``sceptre::import_data`` as ``extra_covariates``. Until this ran, that top-level
frame was empty on the pipeline's own inputs, so SCEPTRE conditioned on nothing
beyond its two automatic depth terms while PerTurbo conditioned on depth, guide
counts and batch.

The library size is deliberately absent: PerTurbo already takes it as an offset and
SCEPTRE adds ``log(response_n_umis)`` itself, so including it again would be
collinear with both.

SCEPTRE drops any factor with fifteen or more levels
(``MAX_N_LEVELS_ALLOWED`` in its formula builder), so on a screen with many
sequencing batches it will silently decline the batch term that PerTurbo uses.
That is a property of SCEPTRE, not something to hide by withholding the column.
"""

from __future__ import annotations

# (column, kind). "continuous" columns go to PerTurbo's --continuous-covariates;
# a "categorical" column goes to --batch-covariate.
CANONICAL_COVARIATES: tuple[tuple[str, str], ...] = (
    ("log1p_total_guide_umis_centered", "continuous"),
    ("percent_mito", "continuous"),
    ("pct_counts_ribo", "continuous"),
    ("batch", "categorical"),
)


def _first_source(mdata, column):
    """The modality frame holding ``column``, preferring the analysed modality."""
    for mod in ("gene", "guide"):
        if mod in mdata.mod and column in mdata[mod].obs.columns:
            return mdata[mod].obs[column]
    return None


def materialize_shared_covariates(mdata, columns=CANONICAL_COVARIATES) -> list[str]:
    """Copy the agreed covariates into the MuData's top-level ``obs``.

    Returns the columns actually written. A column absent from every modality, or
    constant across cells, is skipped: a constant covariate is unidentifiable and
    makes SCEPTRE's regression singular.
    """
    written: list[str] = []
    for column, _kind in columns:
        values = _first_source(mdata, column)
        if values is None:
            continue
        if values.nunique(dropna=False) < 2:
            print(f"Skipping covariate {column!r}: constant across cells.")
            continue
        mdata.obs[column] = values.to_numpy()
        written.append(column)
    print(f"Shared covariates written to the MuData's top-level obs: {written or 'none'}")
    return written


def perturbo_covariate_arguments(available: list[str]) -> tuple[list[str], str | None]:
    """The same covariates, split the way PerTurbo's command line takes them."""
    kinds = dict(CANONICAL_COVARIATES)
    continuous = [c for c in available if kinds.get(c) == "continuous"]
    categorical = [c for c in available if kinds.get(c) == "categorical"]
    return continuous, (categorical[0] if categorical else None)
