# nf-core/crispr: Output

## Introduction

This document describes the output produced by the pipeline.

All paths are relative to the directory supplied with `--outdir`.

The final inference artifacts are written once, under `pipeline_outputs/`. The dashboard directory contains only visualization files and supporting assets; it does not include duplicate copies of the final MuData or local/global analysis TSV outputs.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [Final inference outputs](#final-inference-outputs) - MuData and local/global analysis result tables
- [Pipeline dashboard](#pipeline-dashboard) - Interactive HTML dashboard and plots

### Final inference outputs

<details markdown="1">
<summary>Output files</summary>

- `pipeline_outputs/`
  - `inference_mudata.h5mu`: Final MuData object containing processed modalities and inference results.
  - [`local_analysis_per_element_output.tsv.gz`](https://docs.google.com/spreadsheets/d/1qbXLjDIc5rUJRgl_HhuA68HRxC3Jp5sISr9iqTgyqAM/edit#gid=213615974): Local-analysis element-level inference results.
  - [`local_analysis_per_guide_output.tsv.gz`](https://docs.google.com/spreadsheets/d/1qbXLjDIc5rUJRgl_HhuA68HRxC3Jp5sISr9iqTgyqAM/edit#gid=1700000001): Local-analysis guide-level inference results.
  - [`global_analysis_per_element_output.tsv.gz`](https://docs.google.com/spreadsheets/d/1qbXLjDIc5rUJRgl_HhuA68HRxC3Jp5sISr9iqTgyqAM/edit#gid=864095861): Global-analysis element-level inference results.
  - [`global_analysis_per_guide_output.tsv.gz`](https://docs.google.com/spreadsheets/d/1qbXLjDIc5rUJRgl_HhuA68HRxC3Jp5sISr9iqTgyqAM/edit#gid=1700000002): Global-analysis guide-level inference results.
  - [`catalog_per_element_output.tsv.gz`](https://docs.google.com/spreadsheets/d/1qbXLjDIc5rUJRgl_HhuA68HRxC3Jp5sISr9iqTgyqAM/edit#gid=1841723215): Catalog-formatted per-element table that merges local SCEPTRE and global PerTurbo metrics.
  - [`catalog_per_guide_output.tsv.gz`](https://docs.google.com/spreadsheets/d/1qbXLjDIc5rUJRgl_HhuA68HRxC3Jp5sISr9iqTgyqAM/edit#gid=1700000003): Catalog-formatted per-guide table that merges local SCEPTRE, global PerTurbo, and guide metadata.

</details>

The local-analysis tables contain guide-gene or target-element-gene tests restricted to the configured pairing strategy. The global-analysis tables contain all-by-all PerTurbo tests. Element-level tables aggregate guides by intended target fields: `intended_target_name`, `intended_target_chr`, `intended_target_start`, and `intended_target_end`.

The local/global analysis TSVs include raw p-values, BH-adjusted q-values when applicable, and the nonredundant `*_negLog10p` significance transform. Element-level TSVs also include the catalog-facing element annotations: `element_id`, `element_type`, `element_chr`, `element_start`, `element_end`, `element_name`, `guide_ids`, `num_guides`, `gene_name`, and `nPerturbedCells`. Per-guide TSVs include guide sequence, type, targeting status, guide and intended-target coordinates, PAM/strand, tested-gene symbol, and per-guide perturbed-cell count.

`catalog_per_element_output.tsv.gz` is an additive, per-element catalog view with one row per `(element, gene)` pair and the following columns:
`sceptre_log2_fc`, `sceptre_p_value`, `sceptre_q_value`, `sceptre_fc_se`, `sceptre_negLog10p`, `perturbo_log2_fc`, `perturbo_p_value`, `perturbo_q_value`, `perturbo_fc_se`, `perturbo_negLog10p`,
`element_id`, `element_type`, `element_chr`, `element_start`, `element_end`, `element_name`,
`guide_ids`, `num_guides`, `gene_name`, `gene_id`, and `nPerturbedCells`.

`catalog_per_guide_output.tsv.gz` is the corresponding per-guide catalog view with one row per `(guide_id, gene_id)` pair. It contains the same 10 nonredundant SCEPTRE/PerTurbo metric columns plus:
`guide_id`, `guide_sequence`, `guide_type`, `targeting`, `guide_chr`, `guide_start`, `guide_end`, `guide_strand`, `pam`, `intended_target_name`, `intended_target_chr`, `intended_target_start`, `intended_target_end`, `gene_name`, `gene_id`, and `nPerturbedCells`.

Both catalog schemas are reconstructible from the regular analysis TSVs without opening MuData. For per-element data, outer-join the local SCEPTRE and global PerTurbo columns on `gene_id` plus the intended-target fields. For per-guide data, outer-join them on `(gene_id, guide_id)`; both input tables carry the complete catalog annotation columns.

### Pipeline dashboard

<details markdown="1">
<summary>Output files</summary>

- `pipeline_dashboard/`
  - `dashboard.html`: Interactive dashboard.
  - `figures/`: QC and inference figures used by the dashboard.
  - `evaluation_output/`: Evaluation plots and genome browser files.
  - `guide_seqSpec_plots/`: Guide seqSpec plots.
  - `hashing_seqSpec_plots/`: Hashing seqSpec plots, when hashing is enabled.
  - `additional_qc/`: Additional QC outputs.
  - `benchmark_output/`: Benchmark outputs, when benchmarking is enabled.
- `pipeline_dashboard.tar.gz`: Compressed archive of the `pipeline_dashboard/` directory, written at the top level of `--outdir`.

</details>

`pipeline_dashboard/` intentionally does not contain `inference_mudata.h5mu`, `local_analysis_per_element_output.tsv.gz`, `local_analysis_per_guide_output.tsv.gz`, `global_analysis_per_element_output.tsv.gz`, `global_analysis_per_guide_output.tsv.gz`, `catalog_per_element_output.tsv.gz`, or `catalog_per_guide_output.tsv.gz`. Use the copies in `pipeline_outputs/` as the single source of final analysis outputs.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Resource statistics in `execution_trace.txt`, including requested CPUs/memory/time, attempts, realtime, duration, CPU usage, RSS/VMEM, peak RSS/VMEM, and I/O counters per task.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Original samplesheet resolved from `--input` or the active profile/config: `original_samplesheet.csv` or `original_samplesheet.tsv`.
  - Parameters used by the pipeline run: timestamped `params_*.json` files.
  - Run configuration and log files: `nextflow.config` and `nextflow.log`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
