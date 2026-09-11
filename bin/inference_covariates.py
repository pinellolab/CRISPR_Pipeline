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

The set is deliberately small and follows the runs this method was validated on:
the production Gasperini analysis conditioned on the mitochondrial fraction and
the sequencing batch, and the production Replogle analyses on the batch, both with
the library size as an offset. Applying the same set everywhere matters more than
the exact membership, so both are here and a screen missing one simply skips it.

The library size is absent on purpose: PerTurbo takes it as an offset and SCEPTRE
adds ``log(response_n_umis)`` itself, so listing it again would be collinear with
both. The guide-UMI term the adapter used to pass unilaterally is gone, because
the analyses this method is trusted on did not condition on it and SCEPTRE
excludes gRNA covariates from its own formula by design.

Adding a covariate means adding one line here, and both methods pick it up.

SCEPTRE drops any factor with fifteen or more levels
(``MAX_N_LEVELS_ALLOWED`` in its formula builder), so on a screen with many
sequencing batches it will silently decline the batch term that PerTurbo uses.
That is a property of SCEPTRE, not something to hide by withholding the column.
"""

from __future__ import annotations

# (column, kind, aliases). "continuous" columns go to PerTurbo's --continuous-covariates;
# a "categorical" column goes to --batch-covariate. The aliases are the names the same
# quantity carries in other preprocessing conventions; the first one present is copied
# under the canonical name, so both methods see one column whatever the input called it.
CANONICAL_COVARIATES: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("percent_mito", "continuous", ("percent_mito", "pct_counts_mt", "percent.mito", "pct_mito", "mito_frac")),
    ("batch", "categorical", ("batch",)),
)


def _first_source(mdata, aliases):
    """The first alias present on a modality frame, preferring the analysed modality."""
    for column in aliases:
        for mod in ("gene", "guide"):
            if mod in mdata.mod and column in mdata[mod].obs.columns:
                return column, mdata[mod].obs[column]
    return None, None


def materialize_shared_covariates(mdata, columns=CANONICAL_COVARIATES) -> list[str]:
    """Copy the agreed covariates into the MuData's top-level ``obs``.

    Returns the columns actually written. A column absent from every modality, or
    constant across cells, is skipped: a constant covariate is unidentifiable and
    makes SCEPTRE's regression singular.
    """
    written: list[str] = []
    for column, _kind, aliases in columns:
        source, values = _first_source(mdata, aliases)
        if values is None:
            continue
        if values.nunique(dropna=False) < 2:
            print(f"Skipping covariate {column!r}: constant across cells.")
            continue
        mdata.obs[column] = values.to_numpy()
        # PerTurbo reads the analysed modality's obs by name, so the canonical name
        # has to exist there too when the input used an alias.
        if "gene" in mdata.mod and column not in mdata["gene"].obs.columns:
            mdata["gene"].obs[column] = values.reindex(mdata["gene"].obs_names).to_numpy()
        if source != column:
            print(f"Covariate {column!r} taken from column {source!r}.")
        written.append(column)
    print(f"Shared covariates written to the MuData's top-level obs: {written or 'none'}")
    return written


def perturbo_covariate_arguments(available: list[str]) -> tuple[list[str], str | None]:
    """The same covariates, split the way PerTurbo's command line takes them."""
    kinds = {name: kind for name, kind, _aliases in CANONICAL_COVARIATES}
    continuous = [c for c in available if kinds.get(c) == "continuous"]
    categorical = [c for c in available if kinds.get(c) == "categorical"]
    return continuous, (categorical[0] if categorical else None)
