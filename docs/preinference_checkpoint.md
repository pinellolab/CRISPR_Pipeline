# Pre-inference checkpoint

For a wrapper-free, step-by-step user procedure, see
`docs/preinference_checkpoint_user_guide.md`.

The default `ANALYSIS_mode = 'checkpoint'` runs mapping, RNA/guide/optional HTO
quality control, guide assignment, optional clone removal, and post-concatenation
QC. It then exits successfully **without submitting inference**. The handoff is
published at `OUTDIR/preinference_checkpoint/RUN_NAME/`:

- `preinference_mudata.h5mu`: the unmodified, final pre-inference population;
- `preinference_manifest.json`: its SHA-256, GTF SHA-256, cell/feature counts,
  source commit, and `READY_FOR_INFERENCE` status;
- `checkpoint_qc/`: the published mapping metrics, SeqSpec table, QC tables,
  and plots from the checkpoint. Preserve this directory with the manifest
  and H5MU when moving the handoff.

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

## Dashboard across runs

Open `OUTDIR/pipeline_dashboard/dashboard.html` when the checkpoint completes.
With the execution dashboard enabled, the checkpoint is marked **Ready for
inference**. During continuation, the live execution dashboard shows the
original mapping, preprocessing, guide assignment, and post-concatenation QC
from the copied snapshot under `CONTINUE_OUT/checkpoint_dashboard/`.

After continuation finishes, open
`CONTINUE_OUT/pipeline_dashboard/dashboard.html`. It includes the checkpoint
QC and the resumed run’s inference previews, calibration, and evaluation.
The continuation header records the parent run and the checkpoint, submitted,
and accepted cell counts. Use the accepted count as the input cohort for
inference. The final HTML is also published by standalone Nextflow runs.

For checkpoints created before `checkpoint_qc/` was added, keep the original
output tree in place. The continuation reports whichever published QC plots and
tables it can recover and flags incomplete or unavailable upstream evidence.

`--ANALYSIS_mode full` runs mapping through inference in one invocation.
Edited-H5MU validation currently reads the entire input object into memory.
