# nf-core/crispr: Output

## Introduction

This document describes the output produced by the pipeline.

All paths are relative to the directory supplied with `--outdir`.

## IGVF upload layout (opt-in)

`OUTPUT_IGVF_PACKAGE=true` is the default for a complete `ANALYSIS_mode=full`
run with `INFERENCE_method=default`.
The bundled `nf-crispr-dashboard` plugin is required and is loaded by the pipeline's
default configuration. `OUTPUT_IGVF_TSV=true` (default) streams the four canonical
inference tables from Parquet or compressed TSV to plain TSV inside the pinned base
container. This does not change the inference or its Parquet working outputs.

At successful completion, `results/pipeline_results/` contains exactly the upload
artifacts: `execution_dashboard.html` (the self-contained plugin HTML, including
embedded figures), `inference_mudata.h5mu`, four local/global per-element/per-guide
`*_output.tsv` tables, `nextflow.log`, `samplesheet.csv`,
`full_run_config.json` (resolved parameters, Nextflow version, pipeline Git
commit/branch, command and source manifest), and `pipeline_qc_metrics.json`.
It also embeds the launch configuration sources when their current SHA-256 still
matches the recorded run provenance; credential-setting lines are redacted.
Review any user-supplied parameter values and command-line arguments for
sensitive information before uploading this JSON publicly.
Catalog tables are not part of this upload set. The MuData, log, samplesheet and QC
JSON are hard-linked from their original locations, so they consume no second copy
of the file content. `results/results_extra_files/` contains symlinks to all the
other published top-level result directories/files, preserving existing paths and
avoiding duplication. Upload `pipeline_results/` rather than the symlinked extra
tree. If `OUTPUT_IGVF_TSV=false`, the upload package is incomplete by design;
the original Parquet tables remain available under `pipeline_outputs/`.
Checkpoint runs stop before inference; the current continuation-only and
MuData-only entry paths do not assemble the complete final QC/dashboard bundle,
so they do not claim to create an upload-ready package. Other inference methods
likewise do not have the four canonical local/global tables. These modes retain
their usual outputs without an IGVF package.

For a previously completed run, generate the four TSVs with the pinned base
container, then use `-entry PACKAGE_EXISTING_IGVF_RESULTS --outdir <completed-outdir>`
to create only the package; no mapping or inference runs. On Nextflow 26.04.6,
set `NXF_SYNTAX_PARSER=v1` for the named entry, load the local plugin, and use
`run_support/nextflow_with_storage.sh` so the normal container-storage policy
still applies. The smoke-test config `tests/igvf_packaging_smoke.config` disables
the live observer only for packaging an already-final dashboard; normal runs
keep it enabled.

The final inference artifacts are written once, under `pipeline_outputs/`. The dashboard directory contains visualization files and supporting assets; it does not include duplicate copies of the final MuData or local/global result tables. See the [real result gallery](output-gallery.md) for annotated plots and a completed dashboard package.

See the [MuData field reference](mudata_schema.md) for the internal schema of
`inference_mudata.h5mu` and descriptions of its modalities, matrices,
annotations, assignment layers, and inference-result fields.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps. See [tools and modules](tool-and-module-map.md) for the implementation map and upstream documentation.

| Directory | Main artifact | Use it for |
|---|---|---|
| `pipeline_outputs/` | `inference_mudata.h5mu`, local/global/catalog tables | Reanalysis of cells, guides, and tested pairs. |
| `pipeline_outputs/calibration/` | NTC tables and plots, if enabled and feasible | Checking negative-control behavior. |
| `pipeline_dashboard/` | `execution_dashboard.html`, figures | Browsing workflow status and QC; not the canonical data tables. |
| `pipeline_info/` | Resolved config, provenance, trace, reports | Reproducing the exact run and diagnosing task/resource failures. |

- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution
- [Final inference outputs](#final-inference-outputs) - MuData and local/global analysis result tables
- [Pipeline dashboard](#pipeline-dashboard) - Interactive HTML dashboard and plots

### Final inference outputs

<details markdown="1">
<summary>Output files</summary>

- `pipeline_outputs/`
  - `inference_mudata.h5mu`: Final MuData object containing processed modalities and inference results.
  - `local_analysis_per_element_output.<ext>`: Local-analysis element-level inference results.
  - `local_analysis_per_guide_output.<ext>`: Local-analysis guide-level inference results.
  - `global_analysis_per_element_output.<ext>`: Global-analysis element-level inference results.
  - `global_analysis_per_guide_output.<ext>`: Global-analysis guide-level inference results.
  - `catalog_per_element_output.<ext>`: Catalog-formatted per-element table combining local SCEPTRE and PerTurbo metrics.
  - `catalog_per_guide_output.<ext>`: Catalog-formatted per-guide table combining local SCEPTRE, PerTurbo, and guide metadata.
  - `inference_mode.json`: Resolved inference/control regime.

</details>

`<ext>` is `parquet` by default (`INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT='parquet'`) or `tsv.gz` when that format is selected. The local-analysis tables contain guide-gene or target-element-gene tests restricted to the configured pairing strategy. The global-analysis tables contain all-by-all PerTurbo tests. Element-level tables aggregate guides by intended target fields: `intended_target_name`, `intended_target_chr`, `intended_target_start`, and `intended_target_end`.

The local/global analysis tables include raw p-values, BH-adjusted q-values when applicable, and the nonredundant `*_negLog10p` significance transform. Element-level tables also include the catalog-facing genomic-element annotations: `genomic_element_id`, `genomic_element_type`, `genomic_element_chr`, `genomic_element_start`, `genomic_element_end`, `genomic_element_name`, `guide_ids`, `num_guides`, `gene_name`, and `nPerturbedCells`. Per-guide tables include guide sequence, type, targeting status, guide and intended-target coordinates, PAM/strand, tested-gene symbol, and per-guide perturbed-cell count. Source guide metadata fields named `element_*` are not renamed.

When PerTurbo emits them, its per-pair CRT diagnostics also ride along in the local and global per-guide and per-element result tables under a `perturbo_crt_` prefix (`perturbo_crt_low_information`, `perturbo_crt_observed_nonzero`, `perturbo_crt_expected_nonzero`, `perturbo_crt_saddlepoint_valid`, and the tail-policy columns `perturbo_crt_tail_failure_reason`, `perturbo_crt_used_chernoff`, `perturbo_crt_used_conservative_one`, `perturbo_crt_root_residual_null_sd`). They record how much data a call rests on and how its tail probability was obtained; none of them filters a row, and the catalog tables below do not carry them.

`catalog_per_element_output.<ext>` is an additive, per-element catalog view with one row per `(element, gene)` pair and the following columns:
`sceptre_log2_fc`, `sceptre_p_value`, `sceptre_q_value`, `sceptre_fc_se`, `sceptre_negLog10p`, `perturbo_log2_fc`, `perturbo_p_value`, `perturbo_q_value`, `perturbo_fc_se`, `perturbo_negLog10p`,
`genomic_element_id`, `genomic_element_type`, `genomic_element_chr`, `genomic_element_start`, `genomic_element_end`, `genomic_element_name`,
`guide_ids`, `num_guides`, `gene_name`, `gene_id`, and `nPerturbedCells`.

`catalog_per_guide_output.<ext>` is the corresponding per-guide catalog view with one row per `(guide_id, gene_id)` pair. It contains the same 10 nonredundant SCEPTRE/PerTurbo metric columns plus:
`guide_id`, `guide_sequence`, `guide_type`, `targeting`, `guide_chr`, `guide_start`, `guide_end`, `guide_strand`, `pam`, `intended_target_name`, `intended_target_chr`, `intended_target_start`, `intended_target_end`, `gene_name`, `gene_id`, and `nPerturbedCells`.

Both catalog schemas are reconstructible from the regular analysis tables without opening MuData. For per-element data, outer-join the local SCEPTRE and global PerTurbo columns on `gene_id` plus the intended-target fields. For per-guide data, outer-join them on `(gene_id, guide_id)`; both input tables carry the complete catalog annotation columns.

### Non-targeting calibration and the inference mode

<details markdown="1">
<summary>Output files</summary>

- `pipeline_outputs/inference_mode.json`: the mode both methods received, settled once from the guide assignments: the declared and analysis MOI, the observed assignment counts, SCEPTRE's `control_group`, PerTurbo's `--crt-pool`, and how multi-guide cells are handled.
- `pipeline_outputs/calibration/` (when `ENABLE_NTC_CALIBRATION` is on and both methods ran)
  - `ntc_calibration_per_guide.tsv.gz`, `ntc_calibration_per_element.tsv.gz`: the matched non-targeting pairs with both methods' p-values and log2 fold changes, the pair set, the bucket, cell counts, and PerTurbo's `perturbo_crt_*` diagnostics when present. No q-values: the calibration family is not a discovery family.
  - `ntc_calibration_summary.tsv`, `ntc_calibration_summary.json`: per level, pair set and method, pairs, pairs tested, fractions below 0.01/0.05/0.10, KS distance against the uniform, median p; plus `perturbo_global_all`, PerTurbo's non-targeting rows over every gene.
  - `ntc_method_agreement.tsv`: Pearson on log10 p, Spearman on p, Jaccard at p<0.05 on the pairs both methods tested.
  - `ntc_aucpr_summary.tsv`: direct-target rows against non-targeting rows on the same genes, scored by 1 minus p, AUPRC and AUROC for both methods at both levels, twice: with the negatives drawn 1:1 (`matched_1to1`, the controls evaluation's convention) and with every available negative (`all`), whose precision-recall baseline is the real positive rate. Both methods score the same pairs: a pair either method left without a p-value is dropped for both, and the 1:1 draw is made once and scored by both.
  - `ntc_aucpr_pairs.tsv`: the pairs behind those numbers, per level: every direct-target pair and every non-targeting pair on the scored direct targets' genes, with both p-values, whether the pair was scored and whether it fell in the 1:1 draw.
  - `ntc_cell_count_summary.tsv`, `ntc_calibration_notes.txt`: cells per guide and per element by role, with a warning when the non-targeting median falls outside 0.5x to 2x the targeting median; the regime, the reference pool and the self-inclusion share.
  - Plots: `ntc_calibration_qq_{guide,element}.png`, `ntc_calibration_pvalue_hist_{guide,element}.png`, `ntc_method_scatter_{guide,element}.png`, `ntc_perturbo_global_qq.png`, `ntc_aucpr_curves_{guide,element}.png`, `ntc_cells_per_{guide,element}.png`.
  - `ntc_calibration_skipped.txt` instead, when the screen has no non-targeting guides or, under `nt_cells`, no control-only cell.
- `pipeline_dashboard/benchmark_output/calibration/`: a copy of the directory above.

</details>

### Pipeline dashboard

<details markdown="1">
<summary>Output files</summary>

- `pipeline_dashboard/`
  - `execution_dashboard.html`: Canonical self-contained live/final execution dashboard when the plugin or W&B wrapper is enabled.
  - `dashboard.html`: Compatibility copy of the completed execution dashboard.
  - `legacy_dashboard.html`: Older scientific report retained as an input source, not the supported dashboard.
  - `figures/`: QC and inference figures used by the dashboard.
    - `knee_plot_scRNA_<measurement_set>.png`: independent RNA barcode-rank knee plot for each measurement set.
    - `qc_distributions_scRNA_<measurement_set>.png`: RNA UMI, detected-gene, and mitochondrial distributions with enabled MAD bounds.
    - `rna_qc_filter_flow_<measurement_set>.png`: ordered cells → filter → cells audit diagram with resolved thresholds.
    - `rna_qc_filter_steps_<measurement_set>.png`: before/after metric distributions for every sequential cell filter.
    - `measurement_set_qc_filter_flow.tsv`: exact per-step cell counts, removals, thresholds, and enabled/skipped state.
    - `measurement_set_qc_metrics.tsv`: per-measurement-set cell-retention and QC-threshold audit table.
  - `evaluation_output/`: Evaluation plots and genome browser files.
  - `guide_seqSpec_plots/`: Guide seqSpec plots.
  - `hashing_seqSpec_plots/`: Hashing seqSpec plots, when hashing is enabled.
  - `additional_qc/`: Additional QC outputs.
  - `benchmark_output/`: Benchmark outputs, when benchmarking is enabled.
- `pipeline_dashboard.tar.gz`: Compressed archive of the `pipeline_dashboard/` directory, written at the top level of `--outdir`.

</details>

`pipeline_dashboard/` intentionally does not contain duplicate inference MuData or final local/global/catalog result tables. Use `pipeline_outputs/` as the single source of those analysis outputs.

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
