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
import pandas as pd
import re

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


def declared_external_covariates(mdata) -> list[str]:
    """Return validated numeric covariates introduced at the continuation gate."""
    raw = mdata.uns.get("external_covariates", [])
    if isinstance(raw, (str, bytes)):
        raw = [raw]
    names = [value.decode() if isinstance(value, bytes) else str(value) for value in raw]
    if len(names) != len(set(names)):
        raise ValueError("duplicate external covariate declarations")
    for name in names:
        if not name.startswith("external_") or name not in mdata["gene"].obs:
            raise ValueError(f"invalid or missing external covariate {name!r}")
        values = np.asarray(mdata["gene"].obs[name], dtype=float)
        if not np.isfinite(values).all() or np.unique(values).size < 2:
            raise ValueError(f"external covariate {name!r} must be finite and nonconstant")
    return names


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


def materialize_shared_covariates(mdata, columns=CANONICAL_COVARIATES, config=None) -> list[str]:
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
    if config is not None:
        return _materialize_configured_covariates(mdata, config)
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
    for name in declared_external_covariates(mdata):
        mdata.obs[name] = mdata["gene"].obs[name].reindex(mdata.obs_names).to_numpy()
        written.append(name)
    print(f"Shared covariates written to the MuData's top-level obs: {written or 'none'}")
    return written


def perturbo_covariate_arguments(available: list[str]) -> tuple[list[str], str | None]:
    """The same covariates, split the way PerTurbo's command line takes them."""
    kinds = {name: kind for name, kind, _aliases in CANONICAL_COVARIATES}
    continuous = [c for c in available if kinds.get(c) == "continuous"]
    categorical = [c for c in available if kinds.get(c) == "categorical"]
    return continuous, (categorical[0] if categorical else None)


def _materialize_configured_covariates(mdata, config):
    """Make one full-rank cell design shared by both inference adapters.

    Batch stays categorical; other categorical sample-level terms use reference
    one-hot coding because the pinned PerTurbo CLI accepts only one batch factor.
    """
    gene = mdata["gene"]
    index = gene.obs_names
    candidates = {}
    if config["mito"]:
        _, value = _first_source(mdata, CANONICAL_COVARIATES[0][2])
        if value is None:
            raise ValueError("INFERENCE_COV_MITOCHONDRIAL requested but percent_mito is absent")
        candidates["percent_mito"] = pd.to_numeric(value.reindex(index), errors="raise")
    if config["batch"]:
        _, value = _first_source(mdata, ("batch",))
        if value is None:
            raise ValueError("INFERENCE_COV_BATCH requested but batch is absent")
        candidates["batch"] = value.reindex(index).astype(str)
    if config["ribo"]:
        _, value = _first_source(mdata, ("pct_counts_ribo", "percent_ribo"))
        if value is None:
            raise ValueError("INFERENCE_COV_RIBOSOMAL requested but pre-filter ribosomal percentage is unavailable")
        candidates["percent_ribo"] = pd.to_numeric(value.reindex(index), errors="raise")
    if config["detected_genes"]:
        source, value = _first_source(mdata, ("num_expressed_genes", "n_genes_by_counts", "n_counts", "log1p_n_genes_by_counts"))
        if value is None:
            raise ValueError("INFERENCE_COV_DETECTED_GENES requested but detected gene count is absent")
        measured = pd.to_numeric(value.reindex(index), errors="raise")
        candidates["log1p_detected_genes"] = measured if source == "log1p_n_genes_by_counts" else np.log1p(measured)
    if config["guide_count"]:
        assignment = mdata["guide"].layers.get(ASSIGNMENT_LAYER)
        if assignment is None:
            raise ValueError("INFERENCE_COV_GUIDE_COUNT needs guide_assignment layer")
        candidates["assigned_guide_count"] = pd.Series(
            np.asarray((assignment > 0).sum(axis=1)).ravel(), index=mdata["guide"].obs_names
        ).reindex(index)
    if config["guide_umis"]:
        totals, _ = _guide_umi_totals(mdata)
        candidates["log1p_total_guide_umis"] = pd.Series(
            np.log1p(totals), index=mdata["guide"].obs_names
        ).reindex(index)
    raw_custom = str(config.get("custom", "auto") or "").strip()
    available_custom = sorted(name for name in gene.obs if name.startswith("cov_"))
    if raw_custom.lower() in ("auto", "all"):
        custom = available_custom
    elif raw_custom.lower() in ("", "none"):
        custom = []
    else:
        custom = [name.strip() for name in raw_custom.split(",") if name.strip()]
        missing = sorted(set(custom) - set(available_custom))
        if missing:
            raise ValueError(f"INFERENCE_COV_CUSTOM requested missing columns: {missing}")
    for name in custom:
        values = gene.obs[name].reindex(index)
        if values.isna().any():
            raise ValueError(f"{name} has missing cell values")
        text_values = values.astype(str)
        # IDs with leading zeros are labels, not continuous measurements.
        numeric = pd.to_numeric(text_values, errors="coerce")
        as_numeric = numeric.notna().all() and not text_values.str.match(r"[+-]?0\d+").any()
        if as_numeric:
            candidates[name] = numeric.astype(float)
        else:
            levels = sorted(text_values.unique())
            if len(levels) > 15:
                raise ValueError(f"{name} has {len(levels)} levels; maximum supported is 15")
            encoded_names = set()
            for level in levels[1:]:
                encoded = f"{name}__{re.sub(r'[^A-Za-z0-9_]', '_', level)}"
                if encoded in encoded_names or encoded in candidates:
                    raise ValueError(f"{name} has category labels that collide after safe encoding")
                encoded_names.add(encoded)
                candidates[encoded] = (text_values == level).astype(float)
    for name in declared_external_covariates(mdata):
        candidates[name] = pd.to_numeric(gene.obs[name].reindex(index), errors="raise")

    design = pd.DataFrame(candidates, index=index)
    if design.isna().any().any() or not np.isfinite(design.select_dtypes(include="number").to_numpy()).all():
        raise ValueError("Requested inference covariates contain missing or nonfinite cell values")
    constant = [name for name in design if design[name].nunique(dropna=False) < 2]
    # These two are historical auto-adjustments: single-set demo runs have no
    # estimable batch effect, and a fully MT-filtered tiny pilot can have a
    # constant mitochondrial fraction. Preserve their former skip behavior.
    automatic_skips = [name for name in constant if name in ("batch", "percent_mito")]
    if automatic_skips:
        print(f"Skipping non-estimable automatic covariates: {automatic_skips}")
        design = design.drop(columns=automatic_skips)
    rejected = [name for name in constant if name not in automatic_skips]
    if rejected:
        raise ValueError(f"Requested inference covariates are constant: {rejected}")
    categorical = ["batch"] if "batch" in design else []
    continuous = [name for name in design if name not in categorical]
    rank_parts = [np.ones((len(design), 1), dtype=float)]
    if categorical:
        levels = sorted(design["batch"].astype(str).unique())
        if len(levels) >= 15:
            raise ValueError(f"batch has {len(levels)} levels; SCEPTRE supports fewer than 15")
        rank_parts += [(design["batch"].astype(str) == level).to_numpy(dtype=float)[:, None]
                       for level in levels[1:]]
    rank_parts += [design[name].to_numpy(dtype=float)[:, None] for name in continuous]
    matrix = np.hstack(rank_parts)
    if np.linalg.matrix_rank(matrix) != matrix.shape[1]:
        raise ValueError("Requested covariates are collinear (often cov_* nested in measurement-set batch); disable a redundant term")
    for name in design:
        gene.obs[name] = design[name].to_numpy()
    mdata.obs = design.copy()
    mdata.uns["inference_continuous_covariates"] = continuous
    mdata.uns["inference_batch_covariate"] = categorical[0] if categorical else ""
    mdata.uns["inference_covariate_config"] = [f"{key}={value}" for key, value in sorted(config.items())]
    derive_guide_umi_covariate(mdata)
    print(f"Shared inference covariates: continuous={continuous}, categorical={categorical}")
    return list(design)
