# Configure a run

Start with the basic settings for your assay and execution site. The pipeline loads its root `nextflow.config` automatically. You can pass a small assay config with `-c` to override the relevant settings. In the supplied overlays, active lines set assay values and commented lines show settings you can override.

[Download assay and study configs](examples.md) · [See all parameters](parameters.md)

## Basic settings

| What to check | Parameters and starting values |
|---|---|
| Inputs and results | Set `input` to the [modality-aware samplesheet](samplesheet.md) and `outdir` to a writable results directory. Use paths accessible to every execution node. |
| Run mode | `ANALYSIS_mode = 'checkpoint'` finishes Post Concatenation QC. Use `full` to include inference or `continue` with a saved [checkpoint manifest](preinference_checkpoint_user_guide.md). `DEMO_MODE = true` checks one complete measurement set. |
| Screen design | Set `Multiplicity_of_infection` to `high` or `low`. Review `ENABLE_DATA_HASHING`, `DUAL_GUIDE`, and `is_BaseEditing` against the assay. These three switches default to `false`. |
| RNA counting and barcodes | `scrna_workflow = 'standard'` counts mature transcripts. Use `nac` when the reference includes a nascent index and your assay calls for it. Review `use_multimapping` and `replace_barcodes`, both off by default. Set `bc_replacement_file` when enabling barcode replacement. |
| Guide reads and chemistry | `is_10x3v3` defaults to `false` and controls a feature-barcode mapping override. Review `reverse_complement_guides` and `spacer_tag` against [SeqSpec QC](seqspec.md). A spacer tag enables a whole-read guide search in standard guide mapping. |
| Reference and annotation | The default is `REFERENCE_transcriptome = 'human'` with `use_igvf_reference = true` for the prebuilt RNA index. Use `REFERENCE_gtf_local_path` for an existing annotation or `REFERENCE_gtf_download_path` for the fallback URL. Match the GTF, RNA index, and guide coordinates to the same genome build. |
| Targeted panels | For TAP-seq, use the panel GTF with `REFERENCE_restrict_genes_to_gtf = true`. `TAPSEQ_QC_MODE = true` skips automatic cell-cycle scoring for targeted panels. Review the gene-prevalence cutoff against the genes you expect to measure. |
| RNA and cell QC | Review `QC_barcode_filter` (`knee2`), `QC_min_counts_per_cell` (500 RNA UMIs), `QC_MAD_total_counts` and `QC_MAD_n_genes` (5 each), `QC_pct_mito` (25 percent), and `QC_min_cells_per_gene` (0.05 as a fraction). Each MAD cutoff uses both tails of a log-transformed distribution per measurement set. A value of `0` disables that MAD filter. `ENABLE_SCRUBLET` controls doublet filtering for runs without hashing and defaults to `true`. |
| Guide assignment | `GUIDE_ASSIGNMENT_method = 'sceptre'` is the default. For CLEANSER, set the actual `GUIDE_ASSIGNMENT_capture_method` to `crop-seq` or `direct-capture`. |
| Inference | `INFERENCE_method = 'default'` runs local SCEPTRE and local and global PerTurbo when inference is enabled. `INFERENCE_max_target_distance_bp` defaults to `1000000` (1 Mb) for distance-based local pairs. Global PerTurbo uses its own test universe. |
| Execution site | Select `-profile local`, `slurm`, or `google`. Review `max_cpus`, `max_memory`, the work and container-cache paths, and your scheduler partitions. [Slurm setup](run-recipes.md#slurm-configuration) |

### How knee filtering works

RNA preprocessing automatically creates a barcode-rank plot for each measurement set it processes. `QC_barcode_filter` controls whether the pipeline uses an estimated knee as a cell filter.

| Setting | Barcode filtering |
|---|---|
| `knee2` (default) | Find the first knee, then search the higher-UMI portion of the curve for the second knee. Apply the second knee's RNA UMI threshold when found. |
| `knee` | Apply the first knee's RNA UMI threshold when found. |
| `none` | Skip the knee cutoff. The barcode-rank plot is still generated. |

If the selected knee is unavailable, the pipeline skips that cutoff. The configured RNA UMI minimum and MAD filters still run. CC-Perturb-seq and the Yang example use `none`, while the Gary Hon benchmark example uses `knee`. [QC and filtering](filtering_workflow.md)

## Advanced settings

| When you need more control | Start with |
|---|---|
| Continue a reviewed cohort or run a matched benchmark | `ANALYSIS_parent_manifest`, `ANALYSIS_keep_cells_tsv`, `ANALYSIS_covariates_tsv`, and `BALANCED_*`. [Checkpoint guide](preinference_checkpoint_user_guide.md) |
| Diagnose weak guide recovery or adjust the retained cohort | `GUIDE_MAPPING_QC_*`, assignment probability thresholds, `QC_require_assigned_guide`, guide-count limits, `HTO_*`, and `QC_target_gene_rescue`. [Guide-mapping QC](guide_mapping_overlap_qc.md) · [Target-aware gene filtering](target_aware_gene_filter.md) |
| Change tests or comparison groups | `INFERENCE_target_guide_pairing_strategy`, `INFERENCE_predefined_pairs_to_test`, `INFERENCE_control_group`, and `INFERENCE_COV_*`. [Inference and calibration](inference.md) |
| Tune memory, compute, or model fitting | `INFERENCE_SCEPTRE_*`, `INFERENCE_PERTURBO_*`, and executor or process settings. [All-options Slurm config](downloadable-examples/nextflow_slurm_all_options.config) |
| Change optional diagnostics or output behavior | `QC_EMBEDDING_*`, clone and saturation settings, `NTC_CALIBRATION_*`, `OUTPUT_IGVF_*`, and dashboard settings. [Post Concatenation QC](postconcat_embedding_qc.md) · [Execution dashboard](live_wandb_execution_dashboard.md) |

The [complete parameter reference](parameters.md) lists schema types and defaults. Use the [run recipes](run-recipes.md) for launch commands and record your assay and site overrides with the run.
