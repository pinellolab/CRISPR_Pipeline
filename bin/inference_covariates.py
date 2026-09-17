#!/usr/bin/env python
"""The cell-level covariates both inference methods condition on.

They have to be agreed in one place and written where both methods look, or the
comparison between them is not a comparison of methods. PerTurbo reads named
columns from the analysed modality's ``obs``; SCEPTRE reads the MuData's
top-level ``obs``, which its R reader exposes as ``colData`` and feeds to
``sceptre::import_data`` as ``extra_covariates``.

The columns are written once, in ``bin/mudata_concat.py``, right after the
pipeline's single cell filter and before its gene filter, so every downstream
step (cis subsetting, SCEPTRE chunking, both inference methods) inherits one
set of values. ``ensure_inference_covariates`` is the entry point; it is
idempotent, and ``bin/prepare_inference.py`` and the PerTurbo adapter call it
again only as a fallback for a user-supplied MuData that skipped that step.

The set is the per-cell guide UMI depth, the per-cell library size, the number of
genes detected, and the sequencing batch. All three counts enter on the log
scale.

``percent_mito`` is deliberately absent from this first version.

Why the formula is written out rather than auto-constructed
-----------------------------------------------------------
``bin/inference_sceptre.R`` builds SCEPTRE's formula explicitly instead of
calling ``sceptre:::auto_construct_formula_object``. That builder is a
convenience heuristic with two behaviours that do not suit this pipeline:

* It drops any covariate whose name does not match ``n_umis|n_nonzero`` once the
  column takes fifteen or more distinct values (``MAX_N_LEVELS_ALLOWED``). The
  check has no type test, so it catches continuous covariates as well as
  high-cardinality factors. This is what silently discarded ``percent_mito`` on
  every run, while PerTurbo conditioned on it -- the two methods were never
  fitting the same model.
* It derives the response depth from whatever response matrix SCEPTRE was handed,
  and ``bin/chunk_mudata_sceptre.py`` hands it one gene chunk at a time. On a
  chunked run ``response_n_umis`` is a per-chunk partial sum -- the same cell
  carries a different depth in each chunk -- and ``response_n_nonzero``
  saturates at the chunk size rather than the transcriptome.

SCEPTRE itself has no trouble with continuous covariates: an explicit formula is
passed straight to ``model.matrix``, and a continuous term enters as one numeric
column. The design matrix is still rank-checked --
``convert_covariate_df_to_design_matrix`` raises on a redundant formula -- so
writing the formula out costs no safety.

Writing it out also means these columns can carry their natural names. The
``n_umis`` suffix is only needed to satisfy the auto builder.

One constraint does survive: ``sceptre::import_data`` rejects the names
``response_n_nonzero``, ``response_n_umis``, ``response_p_mito``,
``grna_n_nonzero`` and ``grna_n_umis`` outright, whatever formula is used. The
full-matrix quantities below therefore travel under their own names and replace
SCEPTRE's per-chunk versions in the formula.

SCEPTRE's own gRNA covariates are excluded for a separate reason: at inference
time the gRNA matrix is the binary ``guide_assignment`` layer, so its
``grna_n_umis`` is the number of assigned guides per cell, exactly equal to
``grna_n_nonzero``, and constant at 1 under low MOI. The real depth comes from
``total_guide_umis`` instead.

The guide assignment step needs nothing from here. It runs first
(``guide_assignment_pipeline`` precedes ``inference_pipeline``), it is not
chunked, and no assignment layer exists yet, so SCEPTRE's own
``log(grna_n_umis)`` over raw counts is the genuine depth and
``sceptre::assign_grnas`` already conditions on it. Do not write these columns
into the top-level ``obs`` ahead of that step: ``bin/assign_grnas_sceptre.R``
coerces every ``colData`` column with ``as.factor`` before building a model
matrix, so a continuous column becomes a factor with one level per distinct
value -- a dense ``n x k`` matrix handed to ``Matrix::rankMatrix``.

Two scales, for clarity
-----------------------
SCEPTRE log-transforms inside its formula. PerTurbo standardises every
``--continuous-covariates`` column itself: a column that looks count-like
(non-negative, >= 98% whole numbers) gets ``log1p`` then a z-score, anything
else a z-score alone. Handing it the raw counts would therefore also work, but
which branch fires would depend on a heuristic over the values, and
``log1p`` differs from SCEPTRE's ``log`` at small guide-UMI depths. Writing the
``log`` explicitly makes the scale the same on both sides by construction;
PerTurbo then only centres and scales it, which leaves the fit unchanged. So
each count is written twice: under its natural name for SCEPTRE, and as a
precomputed ``log_`` column for PerTurbo. ``validate_positive`` enforces the
assumption that makes a plain ``log`` safe -- no cell surviving filtering may
have zero depth -- and fails loudly rather than emitting ``-inf``.

What the count columns mean
---------------------------
``total_gene_umis`` is the per-cell UMI total over every gene in the count
matrix and ``num_expressed_genes`` the number of genes with at least one UMI;
both come from scanpy's ``calculate_qc_metrics`` (``total_counts`` and
``n_genes_by_counts``) in ``bin/create_mdata.py``, computed before any gene
restriction. On a targeted panel the matrix *is* the panel, so both are panel
quantities there. ``total_guide_umis`` is the per-cell UMI total over the raw
guide matrix, distinct from ``num_expressed_guides`` (guides with any UMI) and
from the number of *assigned* guides.

``total_gene_umis`` is also PerTurbo's size-factor key (``--library-size-key``),
so PerTurbo receives it both as an offset and as a covariate. That is deliberate:
the offset alone forces a coefficient of one on the log library size, and the
covariate lets the model fit a slope instead.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np

GUIDE_UMI_COLUMN = "total_guide_umis"
GENE_UMI_COLUMN = "total_gene_umis"
DETECTED_GENES_COLUMN = "num_expressed_genes"
BATCH_COLUMN = "batch"
ASSIGNMENT_LAYER = "guide_assignment"

# Names SCEPTRE reserves for the covariates it computes itself. import_data()
# errors if an extra covariate uses one, whatever formula is later supplied.
SCEPTRE_RESERVED_NAMES = frozenset(
    {
        "response_n_nonzero",
        "response_n_umis",
        "response_p_mito",
        "grna_n_nonzero",
        "grna_n_umis",
    }
)


def perturbo_log_name(column: str) -> str:
    return f"log_{column}"


class Covariate(NamedTuple):
    """One conditioned quantity, and how each method reads it.

    ``name`` is the column written to the MuData's top-level ``obs``, which is
    what SCEPTRE sees and what its formula refers to. ``perturbo_name`` is the
    column on the gene modality's ``obs`` that goes on PerTurbo's command line:
    for a count that is the precomputed log, because PerTurbo does not transform.
    """

    name: str
    kind: str  # "count" -> log() for SCEPTRE, log_ column for PerTurbo
    #          # "categorical" -> PerTurbo's --batch-covariate
    aliases: tuple[str, ...]  # input spellings; the first present is used

    @property
    def perturbo_name(self) -> str:
        return perturbo_log_name(self.name) if self.kind == "count" else self.name

    @property
    def sceptre_term(self) -> str:
        """How this covariate is written in SCEPTRE's formula."""
        return f"log({self.name})" if self.kind == "count" else self.name


CANONICAL_COVARIATES: tuple[Covariate, ...] = (
    Covariate(GUIDE_UMI_COLUMN, "count", (GUIDE_UMI_COLUMN,)),
    Covariate(GENE_UMI_COLUMN, "count", (GENE_UMI_COLUMN, "total_counts")),
    Covariate(DETECTED_GENES_COLUMN, "count", (DETECTED_GENES_COLUMN, "n_genes")),
    Covariate(BATCH_COLUMN, "categorical", (BATCH_COLUMN,)),
)

assert not {c.name for c in CANONICAL_COVARIATES} & SCEPTRE_RESERVED_NAMES


def _first_source(mdata, aliases):
    """The first alias present on a modality frame, preferring the analysed modality."""
    for column in aliases:
        for mod in ("gene", "guide"):
            if mod in mdata.mod and column in mdata[mod].obs.columns:
                return column, mdata[mod].obs[column]
    return None, None


def ensure_guide_umi_totals(mdata) -> str | None:
    """Pin ``total_guide_umis`` onto the guide modality while the matrix is whole.

    Call this before any step narrows the guides. The totals are per-cell depth
    over *all* guides, so summing a subset matrix would report a smaller depth --
    and zero for a cell whose only guides were subset away.

    ``bin/create_mdata.py`` normally writes it already. The fallback sums raw
    ``guide.X``, never the ``guide_assignment`` layer: that layer is binary, so
    summing it would give the number of assigned guides rather than UMI depth.
    Nothing in the pipeline binarises ``guide.X`` in place.
    """
    if "guide" not in mdata.mod:
        return None
    guide = mdata["guide"]
    if GUIDE_UMI_COLUMN not in guide.obs.columns:
        guide.obs[GUIDE_UMI_COLUMN] = np.asarray(guide.X.sum(axis=1)).ravel().astype(float)
        print(
            f"Computed {GUIDE_UMI_COLUMN} from guide.X over the full guide matrix "
            f"({guide.n_vars} guides)."
        )
    return GUIDE_UMI_COLUMN


def ensure_gene_depth_totals(mdata) -> list[str]:
    """Pin the cell depth onto the gene modality, before any gene subset.

    ``bin/create_mdata.py`` already writes both columns over every gene in the
    count matrix, so this normally finds them. The fallback sums the matrix it
    is given, so it has to run before the gene filter, the cis subset and the
    SCEPTRE chunking narrow the genes -- summing afterwards would measure depth
    over the retained genes only.
    """
    if "gene" not in mdata.mod:
        return []
    gene = mdata["gene"]
    for column, compute in (
        (GENE_UMI_COLUMN, lambda: np.asarray(gene.X.sum(axis=1)).ravel()),
        (DETECTED_GENES_COLUMN, lambda: np.asarray((gene.X > 0).sum(axis=1)).ravel()),
    ):
        if column not in gene.obs.columns:
            gene.obs[column] = compute().astype(float)
            print(f"Computed {column} over the full gene matrix ({gene.n_vars} genes).")
    return [GENE_UMI_COLUMN, DETECTED_GENES_COLUMN]


def validate_positive(values, column: str) -> np.ndarray:
    """Every count must be strictly positive, so that ``log`` is defined.

    No cell that survives filtering can have zero library size, zero detected
    genes or zero guide UMIs -- ``QC_require_assigned_guide`` in
    ``bin/mudata_concat.py`` alone rules out the last. Rather than paper over a
    violation with ``log1p`` and fit a term that is quietly wrong, fail here and
    say which cells are at fault: it means the filtering upstream did not do what
    this assumes.
    """
    array = np.asarray(values, dtype=float)
    bad = ~(array > 0)
    if bad.any():
        raise ValueError(
            f"Covariate {column!r} must be strictly positive to be log-transformed, "
            f"but {int(bad.sum())} of {array.size} cells are zero, negative or NaN "
            f"(first offending cell index {int(np.flatnonzero(bad)[0])}). "
            "Cells surviving filtering are assumed to have nonzero depth; check the "
            "upstream QC rather than relaxing this."
        )
    return array


def derive_perturbo_log_covariates(mdata) -> list[str]:
    """Write the ``log_`` columns PerTurbo conditions on onto the gene modality.

    PerTurbo applies no transform of its own, so the logs are taken here, over
    the cells of the object it is handed. SCEPTRE gets the untransformed columns
    and writes ``log()`` into its formula instead.
    """
    if "gene" not in mdata.mod:
        return []
    written: list[str] = []
    for covariate in CANONICAL_COVARIATES:
        if covariate.kind != "count":
            continue
        _source, values = _first_source(mdata, covariate.aliases)
        if values is None:
            continue
        array = validate_positive(values.to_numpy(), covariate.name)
        mdata["gene"].obs[covariate.perturbo_name] = np.log(array)
        written.append(covariate.perturbo_name)
    print(f"PerTurbo log covariates written to gene.obs: {written or 'none'}")
    return written


def materialize_shared_covariates(mdata, columns=CANONICAL_COVARIATES) -> list[str]:
    """Copy the agreed covariates into the MuData's top-level ``obs``.

    Returns the SCEPTRE-facing columns actually written. A column absent from
    every modality, or constant across cells, is skipped: a constant covariate is
    unidentifiable and makes the design matrix singular.

    The PerTurbo-scale ``log_`` columns are derived here too, on this object's
    cells, so the analysed MuData carries one set of values rather than the
    adapter's own. They are deliberately not written to the top-level ``obs``:
    handing SCEPTRE both scales of one quantity would be exactly collinear.
    """
    derive_perturbo_log_covariates(mdata)
    written: list[str] = []
    for covariate in columns:
        source, values = _first_source(mdata, covariate.aliases)
        if values is None:
            continue
        if values.nunique(dropna=False) < 2:
            print(f"Skipping covariate {covariate.name!r}: constant across cells.")
            continue
        if covariate.kind == "count":
            validate_positive(values.to_numpy(), covariate.name)
        name = covariate.name
        # Reindex rather than assign positionally: the source may be the guide
        # modality while the destination frame is ordered by the MuData's obs_names.
        mdata.obs[name] = values.reindex(mdata.obs_names).to_numpy()
        if "gene" in mdata.mod and name not in mdata["gene"].obs.columns:
            mdata["gene"].obs[name] = values.reindex(mdata["gene"].obs_names).to_numpy()
        if source != name:
            print(f"Covariate {name!r} taken from column {source!r}.")
        written.append(name)
    print(f"Shared covariates written to the MuData's top-level obs: {written or 'none'}")
    return written


def ensure_inference_covariates(mdata) -> list[str]:
    """Derive and write every conditioned covariate, once, on the object given.

    The totals are pinned first (they must see the whole matrices), then the
    ``log_`` columns for PerTurbo and the plain columns for SCEPTRE. Calling it
    again on an object that already carries the columns rewrites the same
    values, so the fallbacks downstream are safe.
    """
    ensure_guide_umi_totals(mdata)
    ensure_gene_depth_totals(mdata)
    return materialize_shared_covariates(mdata)


def sceptre_formula_terms(available: list[str]) -> list[str]:
    """The formula terms SCEPTRE should fit, for the covariates actually present.

    ``bin/inference_sceptre.R`` builds the same list on its own side; this exists
    so the Python tests can assert the two agree.
    """
    return [c.sceptre_term for c in CANONICAL_COVARIATES if c.name in set(available)]


def perturbo_covariate_arguments(available) -> tuple[list[str], str | None]:
    """The same covariates, split the way PerTurbo's command line takes them.

    The names returned are the PerTurbo-side columns: the precomputed logs for
    the counts, the plain column for the batch.
    """
    available = list(available)
    continuous = [c.perturbo_name for c in available if c.kind == "count"]
    categorical = [c.perturbo_name for c in available if c.kind == "categorical"]
    return continuous, (categorical[0] if categorical else None)
