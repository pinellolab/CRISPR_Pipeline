# Pause before inference, edit cells or covariates, then continue

This guide applies to the current `dev` pipeline. It uses
**two ordinary Nextflow invocations**, not a custom run wrapper. The first
ends successfully at a pre-inference MuData checkpoint; the second starts
inference from that checkpoint. `checkpoint` is the default `ANALYSIS_mode`; explicitly select `full` for unattended inference.

## The distinction that matters

- `--ANALYSIS_mode checkpoint` is the intentional **stop point**. Mapping,
  preprocessing, guide assignment, and configured pre-inference QC finish;
  inference is not scheduled. The run status is successful, not paused or
  failed. This branch defaults to `checkpoint`, but set it explicitly.
- `--ANALYSIS_mode continue` is a **new Nextflow run**. It validates a
  checkpoint and optional external edits, then runs SCEPTRE, PerTurbo, and
  evaluation. It does not remap reads or redo guide assignment.
- Nextflow's `-resume` has a different purpose: it retries an interrupted
  invocation using that invocation's history and work directory. It does not
  turn a checkpoint run into an inference run.
- Use a different output directory for every continuation. Keep the parent
  checkpoint and all external input files immutable.

## 1. Prerequisites

Have a clean checkout of the pinned pipeline commit, Nextflow, the dataset's existing
params file, its samplesheet with accessible local paths, any resource config,
and the same reference GTF for both stages. Pin container images in the
dataset config. Review assay-specific parameters before the first launch.
The first pilot supports default inference method and default
target–guide pairing on continuation.

The following variables only shorten the displayed commands; replace their
values with your own paths and run names. Use a persistent, shared work
directory if your executor uses multiple nodes.

```bash
PIPELINE=/absolute/path/CRISPR_Pipeline
PARAMS=/absolute/path/dataset.params.json
SAMPLESHEET=/absolute/path/samplesheet.local.csv
RESOURCES=/absolute/path/dataset.resources.config
GUIDES=/absolute/path/guide_metadata.tsv
GTF=/absolute/path/reference.gtf.gz
PARENT_OUT=/absolute/path/results_before_external_edit
CONTINUE_OUT=/absolute/path/results_after_external_edit
PARENT_WORK=/absolute/path/work/checkpoint
CONTINUE_WORK=/absolute/path/work/continue
PARENT_NAME=sample_checkpoint_01
CONTINUE_NAME=sample_continue_01
```

The checkout's `nextflow.config` supplies the remaining defaults. Do not
change the checkout between a run's provenance preflight and its launch.

## 2. Run to the checkpoint

From a stable launch directory, run:

```bash
nextflow -c "$RESOURCES" run "$PIPELINE/main.nf" \
  -profile local -params-file "$PARAMS" \
  --input "$SAMPLESHEET" --outdir "$PARENT_OUT" \
  --REFERENCE_gtf_local_path "$GTF" \
  --ANALYSIS_mode checkpoint \
  -work-dir "$PARENT_WORK" -name "$PARENT_NAME"
```

For a cluster, replace `local` and the resource config with the executor
profile and shared paths appropriate to that cluster. The result to preserve
is:

```text
$PARENT_OUT/preinference_checkpoint/$PARENT_NAME/
  preinference_mudata.h5mu
  preinference_manifest.json
```

Check that the run completed successfully and the manifest says
`READY_FOR_INFERENCE`. The manifest records parent cell/feature counts and
hashes of the MuData and reference GTF. **Do not edit the parent H5MU or
manifest.**

If the checkpoint invocation fails before publishing the handoff, retry that
same stage from its original launch directory with the same work directory
and `-resume`. Some Nextflow installations reject reusing a completed run
name; choose a new attempt name if necessary and use the manifest in the
actual published run-name directory. Do not use `continue` until the
checkpoint exists and is ready.

## 3. Prepare the external selection

The simplest handoff needs **no edited H5MU**. Create a tab-separated
allowlist with exactly one header, `cell_id`, containing the *full, exact*
cell identifiers in the parent MuData:

```text
cell_id
CELL_A-1_IGVFDSXXXX
CELL_B-1_IGVFDSXXXX
```

If adding covariates, create another tab-separated file keyed by `cell_id`:

```text
cell_id	external_treatment_score
CELL_A-1_IGVFDSXXXX	0.25
CELL_B-1_IGVFDSXXXX	-0.40
```

The cell IDs above are **illustrative**, not IDs from a real dataset. Use
the identifiers in `preinference_mudata.h5mu`, including measurement-set
suffixes. Each covariate name must start with `external_`; pilot support is
limited to numeric, finite, nonconstant values. Supply a value for every
retained cell. Do not add a covariate merely because it is available: choose
one justified by the experiment and avoid redundant or outcome-derived terms.

If an external tool must edit MuData itself, write a **new** H5MU and pass it
as `--ANALYSIS_input_mudata`. It may remove cells but may not add cells or
change gene/guide feature IDs, raw gene/guide counts, or the guide-assignment
layer. The continuation validates these constraints. You can supply an
allowlist and covariate TSV with this edited H5MU as well; the allowlist must
refer only to cells remaining in the edited file.

Put each proposed edit in a new, versioned directory. Do not overwrite a TSV
or H5MU used by an earlier continuation: that makes its provenance ambiguous
and can make Nextflow cache reuse hard to reason about.

## 4. Continue into inference

With an allowlist and numeric covariates, but **without** editing the H5MU:

```bash
MANIFEST="$PARENT_OUT/preinference_checkpoint/$PARENT_NAME/preinference_manifest.json"
KEEP=/absolute/path/edits_v1/keep_cells.tsv
COVARIATES=/absolute/path/edits_v1/covariates.tsv

nextflow -c "$RESOURCES" run "$PIPELINE/main.nf" \
  -profile local -params-file "$PARAMS" \
  --outdir "$CONTINUE_OUT" --ANALYSIS_mode continue \
  --REFERENCE_gtf_local_path "$GTF" \
  --ANALYSIS_parent_manifest "$MANIFEST" \
  --ANALYSIS_keep_cells_tsv "$KEEP" \
  --ANALYSIS_covariates_tsv "$COVARIATES" \
  -work-dir "$CONTINUE_WORK" -name "$CONTINUE_NAME"
```

If you also created an edited H5MU, add
`--ANALYSIS_input_mudata /absolute/path/edits_v1/edited.h5mu`. If you only
want filtering, omit `--ANALYSIS_covariates_tsv`; if you only want covariates,
omit `--ANALYSIS_keep_cells_tsv`. Omit all three optional edit arguments to
continue with the original checkpoint unchanged.

Use the **same reference GTF content** as the checkpoint stage. The
continuation compares its SHA-256 against the parent manifest and refuses a
mismatch. A params file may carry the same local GTF path for both stages;
there is no need to remap FASTQs in the continuation.

The continuation publishes:

```text
$CONTINUE_OUT/external_mudata/analysis_input.h5mu
$CONTINUE_OUT/external_mudata/external_mudata_report.json
$CONTINUE_OUT/pipeline_outputs/inference_mudata.h5mu
$CONTINUE_OUT/pipeline_outputs/*analysis*_output.parquet
```

Read `external_mudata_report.json` before interpreting inference. It records
parent, submitted, and accepted cell counts, external covariate names, and
file hashes. Check the Nextflow task list: this run should contain
`prepare_external_mudata`, inference, and evaluation, but no mapping task.

If **this** invocation fails after starting inference, retry from its own
launch directory with its own work directory and `-resume`; keep the parent
manifest and edit files identical. Changing the cell set or covariates is a
**new continuation**, not a retry: use new immutable edit files, output
directory, work directory, and run name.

## Local Pinello-server preflight (required on this server)

The two portable commands above omit host-specific storage and provenance.
For **each** new or resumed local run on `/data/pinello/lucas_logan`, keep
the checkout clean and do the preflight inline—no wrapper script is needed.
For the checkpoint stage:

```bash
PARENT_PROV="$PARENT_OUT/pipeline_info/$PARENT_NAME"
python "$PIPELINE/bin/run_provenance_agent.py" prepare \
  --pipeline-repo "$PIPELINE" --source-dir "$PIPELINE" \
  --output-dir "$PARENT_PROV" --run-name "$PARENT_NAME" \
  --artifact "params=$PARAMS" --artifact "resources=$RESOURCES" \
  --artifact "samplesheet=$SAMPLESHEET" --artifact "guide_metadata=$GUIDES" \
  --artifact "gtf=$GTF" \
  --artifact "storage_helper=/data/pinello/lucas_logan/run_support/container_storage.sh" \
  --artifact "storage_config=/data/pinello/lucas_logan/run_support/container_storage.config"
python "$PIPELINE/bin/run_provenance_agent.py" check \
  --pipeline-repo "$PIPELINE" \
  --metadata "$PARENT_PROV/repository_provenance.json"
```

For the continuation stage, after creating `MANIFEST`, `KEEP`, and
`COVARIATES`:

```bash
CONTINUE_PROV="$CONTINUE_OUT/pipeline_info/$CONTINUE_NAME"
python "$PIPELINE/bin/run_provenance_agent.py" prepare \
  --pipeline-repo "$PIPELINE" --source-dir "$PIPELINE" \
  --output-dir "$CONTINUE_PROV" --run-name "$CONTINUE_NAME" \
  --artifact "params=$PARAMS" --artifact "resources=$RESOURCES" \
  --artifact "guide_metadata=$GUIDES" --artifact "gtf=$GTF" \
  --artifact "parent_manifest=$MANIFEST" \
  --artifact "parent_mudata=$PARENT_OUT/preinference_checkpoint/$PARENT_NAME/preinference_mudata.h5mu" \
  --artifact "cell_allowlist=$KEEP" --artifact "covariates=$COVARIATES" \
  --artifact "storage_helper=/data/pinello/lucas_logan/run_support/container_storage.sh" \
  --artifact "storage_config=/data/pinello/lucas_logan/run_support/container_storage.config"
python "$PIPELINE/bin/run_provenance_agent.py" check \
  --pipeline-repo "$PIPELINE" \
  --metadata "$CONTINUE_PROV/repository_provenance.json"
```

Add any other dataset-specific files (such as a metadata repair report or
edited H5MU) as named `--artifact` arguments; omit artifacts for optional
files you did not create. For the actual local launch, insert these two
additional config arguments **after** `-c "$RESOURCES"` and before `run`
in the corresponding command above:

```bash
-c /data/pinello/lucas_logan/run_support/container_storage.config \
-c "$PARENT_PROV/provenance.generated.config"
```

Use `$CONTINUE_PROV` instead of `$PARENT_PROV` for the continuation. Source
`/data/pinello/lucas_logan/run_support/container_storage.sh` **immediately
before** each Nextflow launch and stop if it returns nonzero. Do not redirect
image cache or unpack paths to home or system `/tmp` afterward. Preserve the
generated provenance files with each run. This policy is specific to this
server; remote SLURM/cloud nodes need their own accessible storage paths.

## Current limitations

- The continuation runs inference/evaluation but does **not**
  reconstruct the parent mapping/QC HTML dashboard or publish a combined
  before/after cell-flow card to W&B. The external-MuData report is the
  authoritative count comparison for now.
- The edited-H5MU validator loads both MuData objects into memory; for very
  large experiments, use an allowlist/covariate TSV when possible and plan
  sufficient memory for validation.
- Only numeric external covariates are supported by this continuation interface.
- This is not an automatic human-in-the-loop pause within one Nextflow
  session. It is a deliberate successful checkpoint run followed by a
  separately launched continuation.
