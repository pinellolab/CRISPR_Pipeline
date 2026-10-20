# Pre-inference checkpoint pilot

The default `ANALYSIS_mode = 'checkpoint'` runs mapping, RNA/guide/optional HTO
quality control, guide assignment, optional clone removal, and post-concatenation
QC. It then exits successfully **without submitting inference**. The handoff is
published at `OUTDIR/preinference_checkpoint/`:

- `preinference_mudata.h5mu`: the unmodified, final pre-inference population;
- `preinference_manifest.json`: its SHA-256, GTF SHA-256, cell/feature counts,
  source commit, and `READY_FOR_INFERENCE` status.

Run a continuation with a *different output directory* and
`--ANALYSIS_mode continue --ANALYSIS_parent_manifest /path/to/preinference_manifest.json`.
Without optional edits, it uses the parent MuData unchanged. Optional edits are:

- `--ANALYSIS_input_mudata /path/to/edited.h5mu`: a copied H5MU that may contain
  fewer cells, but cannot add cells or alter RNA counts, guides, or assignments;
- `--ANALYSIS_keep_cells_tsv /path/to/keep.tsv`: a tab-separated, one-column
  `cell_id` allowlist, matching the canonical MuData cell IDs exactly;
- `--ANALYSIS_covariates_tsv /path/to/covariates.tsv`: a tab-separated table
  keyed by `cell_id`, with one or more finite, nonconstant **numeric** columns
  named `external_*`. This first pilot rejects categorical custom covariates.

The continuation validates the parent and reference hashes, applies the
external selection, and publishes `external_mudata/analysis_input.h5mu` and
`external_mudata/external_mudata_report.json` before invoking the existing
SCEPTRE/PerTurbo inference workflow. Both methods receive each declared numeric
external covariate. Evaluation runs afterward. The original checkpoint is
never overwritten. `-resume` retries tasks within a given invocation; changing
an edited H5MU in place is not a supported cache-invalidation mechanism.

`--ANALYSIS_mode full` explicitly opts into the previous uninterrupted path.
This remains a feature-branch pilot: the continuation does not yet reconstruct
the full mapping/assignment HTML dashboard from parent-run artifacts, and the
edited-H5MU validation currently reads the entire object into memory.
