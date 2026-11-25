# Run recipes

These commands show the supported entry points. Replace illustrative paths and select the executor/profile configured for your site. Keep a stable launch directory, work directory, source commit, params/configs, and inputs for every `-resume`. On the Pinello server, also perform the mandatory [container-storage and provenance preflight](preinference_checkpoint_user_guide.md#local-pinello-server-preflight-required-on-this-server) immediately before **each** launch; do not use home or system `/tmp` as image caches.

## Slurm configuration

Use the immutable `v1.1.0` tag for a reproducible deployment. The
[all-options Slurm template](downloadable-examples/nextflow_slurm_all_options.config)
lists active defaults and site settings by topic. The
[compact Yang overlay](downloadable-examples/nextflow_yang.dev.config)
contains the Yang assay and QC overrides and defaults to a full analysis.
Choose one overlay and supply your own samplesheet, output directory, shared
work and cache paths, and scheduler partitions.

```bash
git clone --branch v1.1.0 https://github.com/pinellolab/CRISPR_Pipeline.git
cd CRISPR_Pipeline
export NXF_VER=26.04.6
export NXF_PLUGINS_DIR=/shared/nextflow-plugins
export NXF_SYNTAX_PARSER=v1
bash run_support/install_dashboard_plugin.sh
nextflow -version
export YANG_SLURM_QUEUE=your_cpu_partition
export YANG_SIF_CACHE=/shared/singularity-cache
nextflow -c conf/simplified_pipeline/nextflow_yang.dev.config run main.nf \
  -profile slurm -work-dir /shared/work/yang \
  --input /shared/inputs/samplesheet.csv --outdir /shared/results/yang_qc \
  --ANALYSIS_mode checkpoint
```

Run these exports in the environment that starts Nextflow, including inside
your launcher `sbatch` script if you use one. For a non-Yang study, use the
[matching assay overlay](examples.md) with a separate site config. For example,
set the CPU queue using `profiles { slurm { process.queue = 'your_cpu_partition' } }`
in a site config and supply it as another `-c` file.

Use `checkpoint` to finish dashboard Step 7 and write the post-QC MuData.
For a GPU-enabled full run, set the optional GPU partition and any required
site-specific Slurm GRES value in your copy of the overlay. See
[plugin installation](plugin-installation.md) for local image paths and
dashboard configuration.

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

## Google Batch

The `google` profile uses Google Batch. Set `google_project`, `google_region`
and `google_bucket` for your deployment, and authenticate the Nextflow launcher
with credentials for your project. The Batch workers need access to the input
objects, container images and work bucket. Use the same [dataset config](examples.md)
and run mode with `-profile google`. Use site-specific container and
compute overrides for workloads such as GPU PerTurbo. [Nextflow Google Batch documentation](https://docs.seqera.io/nextflow/google)

## Validate and size your configuration

Use `nextflow -c /path/site.config config . -profile slurm -flat` to review
the effective executor, queue, container settings and assay parameters. The
same command with `-profile local` checks a workstation configuration.

Treat `max_cpus` and `max_memory` as per-task limits used by resource
selectors. A Slurm partition must accommodate each task's actual request.
Use a small measurement set and inspect the trace before setting limits for
a full screen. Large RNA mapping, concatenation and global inference often
have different bottlenecks. For run-specific task settings and work
directories, see [outputs and provenance](output.md).

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

All options and defaults are in the [configuration guide](configuration.md). Record scientific and computational choices in the run's resolved configuration.
