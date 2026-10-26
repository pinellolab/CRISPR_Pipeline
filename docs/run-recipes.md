# Run recipes

These commands show the supported entry points. Replace illustrative paths and select the executor/profile configured for your site. Keep a stable launch directory, work directory, source commit, params/configs, and inputs for every `-resume`. On the Pinello server, also perform the mandatory [container-storage and provenance preflight](preinference_checkpoint_user_guide.md#local-pinello-server-preflight-required-on-this-server) immediately before **each** launch; do not use home or system `/tmp` as image caches.

## End-to-end from FASTQs

The default is a pre-inference checkpoint, so set `full` explicitly when you want inference in the same invocation:

```bash
nextflow run /path/CRISPR_Pipeline/main.nf \
  -profile local -c /path/dataset.resources.config \
  -params-file /path/dataset.params.json \
  --input /path/samplesheet.csv --outdir /path/results_full \
  --ANALYSIS_mode full \
  -work-dir /path/work_full -name dataset_full_01
```

Use a pinned, compatible PerTurbo image from the current pipeline config. `INFERENCE_method=default` runs local SCEPTRE plus local/global PerTurbo. The global result format defaults to Parquet, which needs the pinned image's Parquet support. Check the resolved config and container provenance in `pipeline_info`.

## A small mapping/QC preflight

Add `--DEMO_MODE true` to select the eligible measurement set with the fewest RNA FASTQ files; all lanes and modalities for that set are retained. Use a separate output/work directory and a clearly marked run name. This is not a final biological analysis.

## Checkpoint, edit cells or covariates, continue

Run the first command with `--ANALYSIS_mode checkpoint`. It publishes `preinference_checkpoint/RUN_NAME/preinference_mudata.h5mu` and `preinference_manifest.json` without inference. Then create immutable `keep_cells.tsv` and/or `covariates.tsv` files—or a copied edited H5MU—and start a **new** invocation with `--ANALYSIS_mode continue --ANALYSIS_parent_manifest /path/to/preinference_manifest.json`. The [step-by-step checkpoint guide](preinference_checkpoint_user_guide.md) includes exact commands, validation rules, output locations, and the distinction from `-resume`.

## Retry a failed invocation

Use the same launch directory, source commit, input files, params, work directory and Nextflow session; add `-resume`. If an input changed, expect affected tasks to rerun. Never overwrite a checkpoint or external edit file in place to trick the cache. `-resume` does not skip QC or turn a checkpoint run into a continuation.

## Start inference from an existing MuData

The separate `INFERENCE_FROM_MUDATA` entrypoint accepts a post-guide-assignment H5MU and reruns default inference, but does **not** validate its lineage against a checkpoint manifest or reconstruct mapping/QC dashboards. Prefer `ANALYSIS_mode=continue` for auditable external filtering.

```bash
nextflow run /path/CRISPR_Pipeline/main.nf \
  -entry INFERENCE_FROM_MUDATA -profile local \
  -c /path/dataset.resources.config \
  --INFERENCE_input_mudata /path/input.h5mu \
  --INFERENCE_method default \
  --INFERENCE_target_guide_pairing_strategy default \
  --REFERENCE_gtf_local_path /path/reference.gtf.gz \
  --outdir /path/inference_only -work-dir /path/work_inference_only
```

## Technology switches

Put assay-specific settings in a reviewed dataset config or params file rather than editing `nextflow.config`:

| Assay | Key choices to verify |
|---|---|
| HTO/multiplexed | `ENABLE_DATA_HASHING=true`, `hash` rows, HTO metadata and SeqSpec, barcode replacement if the modalities use different barcode representations. HTO singlet filtering replaces Scrublet in this path. |
| TAP-seq targeted panel | `TAPSEQ_QC_MODE=true`, an appropriate panel GTF, `REFERENCE_restrict_genes_to_gtf=true`, and a biologically justified small `QC_min_cells_per_gene` fraction. |
| Dual-guide | `DUAL_GUIDE=true`, guide pairs and non-targeting controls represented in metadata; inspect the post-collapse cell count. |
| Base editing | `is_BaseEditing=true` and, for the streaming mapper, `BASEEDITING_method='flash'`; inspect guide length/tolerance/orientation against reads. |
| Low MOI | `Multiplicity_of_infection='low'`; verify non-targeting controls exist for the default control-anchored contrast. |
| Guide calling | `GUIDE_ASSIGNMENT_method='sceptre'` or `'cleanser'`; for CLEANSER choose the actual capture method (`crop-seq` or `direct-capture`). |

All options and defaults are in the [generated parameter reference](parameters.md). Scientific and computational choices should be recorded in the run's resolved configuration; run-specific parameters are not interchangeable across assays.
