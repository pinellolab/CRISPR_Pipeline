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

One derived column is defined here without being in that list:
``log1p_total_guide_umis_centered`` (see ``derive_guide_umi_covariate``). It is
not a member of ``CANONICAL_COVARIATES``, so neither method conditions on it --
that is the decision recorded in the paragraph above and it is unchanged. What
is defined here is the *derivation*, in one function, because the column is
still written into the analysed object and read back by the PerTurbo adapter's
zero-variance guard, and two copies of the formula centred over two cell
populations is exactly the kind of difference this module exists to remove.
The centring mean is taken over the cells in the object it is given, which is
now the same population for both methods (``QC_require_assigned_guide`` in
``bin/mudata_concat.py``).

SCEPTRE drops any factor with fifteen or more levels
(``MAX_N_LEVELS_ALLOWED`` in its formula builder), so on a screen with many
sequencing batches it will silently decline the batch term that PerTurbo uses.
That is a property of SCEPTRE, not something to hide by withholding the column.
"""

from __future__ import annotations

import numpy as np

GUIDE_UMI_COLUMN = "total_guide_umis"
GUIDE_UMI_COVARIATE = "log1p_total_guide_umis_centered"
ASSIGNMENT_LAYER = "guide_assignment"

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


def _guide_umi_totals(mdata):
    """Per-cell guide UMI totals, from the column the pipeline already computes.

    ``bin/create_mdata.py`` writes ``guide.obs['total_guide_umis']`` from the full
    guide matrix, so it does not move when a later step subsets the guides. Only
    when it is absent are the totals summed here, preferring the assignment layer
    over raw ``guide.X`` -- the same order the PerTurbo adapter used.
    """
    guide = mdata["guide"]
    if GUIDE_UMI_COLUMN in guide.obs.columns:
        return np.asarray(guide.obs[GUIDE_UMI_COLUMN], dtype=float), GUIDE_UMI_COLUMN
    matrix = guide.layers[ASSIGNMENT_LAYER] if ASSIGNMENT_LAYER in guide.layers else guide.X
    source = (
        f"guide.layers['{ASSIGNMENT_LAYER}']" if ASSIGNMENT_LAYER in guide.layers else "guide.X"
    )
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    guide.obs[GUIDE_UMI_COLUMN] = totals
    return totals, source


def derive_guide_umi_covariate(mdata) -> str | None:
    """Write ``log1p_total_guide_umis_centered`` onto the gene modality's ``obs``.

    One definition of the derivation, centred over the cells of the object it is
    handed. Neither method conditions on it -- it is not in
    ``CANONICAL_COVARIATES`` -- but it is written into the analysed object and the
    PerTurbo adapter reads it back for its zero-variance-in-controls guard, so it
    must not be computed twice over two different cell populations.

    Returns the column name, or ``None`` when there is no guide modality to
    derive it from.
    """
    if "guide" not in mdata.mod or "gene" not in mdata.mod:
        return None
    totals, source = _guide_umi_totals(mdata)
    centered = np.log1p(totals)
    centered = centered - float(np.nanmean(centered))
    mdata["gene"].obs[GUIDE_UMI_COVARIATE] = centered
    print(
        f"Derived {GUIDE_UMI_COVARIATE} from {source} over "
        f"{mdata['gene'].n_obs} cells (not a conditioned covariate)."
    )
    return GUIDE_UMI_COVARIATE


def materialize_shared_covariates(mdata, columns=CANONICAL_COVARIATES) -> list[str]:
    """Copy the agreed covariates into the MuData's top-level ``obs``.

    Returns the columns actually written. A column absent from every modality, or
    constant across cells, is skipped: a constant covariate is unidentifiable and
    makes SCEPTRE's regression singular.

    The derived guide-UMI column is computed here too, on this object's cells, so
    that the analysed MuData carries one set of values rather than the adapter's
    own. It is deliberately not written to the top-level ``obs`` and not returned:
    it is not in ``CANONICAL_COVARIATES``, so neither method conditions on it, and
    putting it in SCEPTRE's ``extra_covariates`` would silently add a gRNA
    covariate SCEPTRE excludes by design.
    """
    derive_guide_umi_covariate(mdata)
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
