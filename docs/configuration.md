# Configure a run

[Basic settings](#basic-settings) · [Config downloads](#config-downloads) · [Advanced settings](#advanced-settings)

The pipeline loads its root `nextflow.config` automatically. Add a study or site overlay with `-c` and set run-specific values on the command line, in the overlay, or in a JSON or YAML params file. The tables show configured defaults on `main`. A tagged checkout or an overlay can supply different values. [Run recipes](run-recipes.md)

## Basic settings

For a FASTQ run, pass `--input /path/to/samplesheet.csv` and `--outdir /path/to/results` at launch. Use a [modality-aware samplesheet](samplesheet.md) and paths accessible to every compute node. Pick an executor with `-profile local`, `-profile slurm`, or `-profile google`. [Slurm setup](run-recipes.md#slurm-configuration)

<!-- BASIC_PARAMETERS -->

### How knee filtering works

RNA preprocessing creates a barcode-rank plot for each measurement set. `knee`, the default on `main`, applies the first estimated knee. `knee2` searches for a second knee in the higher-UMI portion of the curve. `none` keeps the plot and skips the knee cutoff. If the selected knee is unavailable, the configured RNA UMI minimum and MAD filters still run. [QC and filtering](filtering_workflow.md)

The `v1.1.0` tag defaults to `knee2`. The Gasperini, Replogle essential-gene, and TAP-seq examples explicitly use `knee2`. The Gary Hon benchmark example uses `knee`. CC-Perturb-seq and the Yang example use `none`.

## Config downloads

In an overlay, active lines specify run or assay values and commented lines show settings you can override. Keep the root config in its pipeline checkout so it can resolve its included files.

| Config | Use |
|---|---|
| [nextflow.config](https://github.com/pinellolab/CRISPR_Pipeline/blob/main/nextflow.config) | Root config with run defaults and executor profiles. |
| [nextflow_cc.config](downloadable-examples/nextflow_cc.config) | CC-Perturb-seq with barcode replacement. |
| [nextflow_tapseq.config](downloadable-examples/nextflow_tapseq.config) | Targeted TAP-seq. Supply the panel GTF via `--REFERENCE_gtf_local_path`. |
| [nextflow_yang.dev.config](downloadable-examples/nextflow_yang.dev.config) | Compact Yang assay and Slurm overlay. |
| [nextflow_slurm_all_options.config](downloadable-examples/nextflow_slurm_all_options.config) | Commented Slurm template with basic and advanced options and site settings. |
| [Gasperini config](downloadable-examples/gasperini_at_scale/nextflow.config) | Gasperini at-scale K562 CRISPRi. |
| [Replogle config](downloadable-examples/replogle_k562_essential/nextflow.config) | Replogle K562 day-6 essential genes. |
| [Gary Hon config](downloadable-examples/gary_hon_tf_benchmark/nextflow.config) | Gary Hon WTC11 TF benchmark with HTOs. |

The study bundles include the samplesheets, metadata, and SeqSpecs used with these configs. [Study inputs](examples.md#study-inputs)

## Advanced settings

For config-only launches, set `params.input` in your run config. The base config leaves it `null`. Continuations use the saved checkpoint manifest. [Checkpoint guide](preinference_checkpoint_user_guide.md)

Use the sections below to find when each setting applies. [Inference and calibration](inference.md) · [Post Concatenation QC](postconcat_embedding_qc.md) · [Execution dashboard](live_wandb_execution_dashboard.md)

<!-- ADVANCED_PARAMETERS -->
