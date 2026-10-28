# Inference, covariates, and calibration

The default inference runs local/cis SCEPTRE, local/cis PerTurbo and global/trans PerTurbo from the same final raw-count, guide-assigned MuData. The current PerTurbo local path passes explicit `--pairs-to-test`; the global path tests its full configured universe. A pair absent from local results may have been excluded by distance, missing guide/element metadata, gene filtering, or insufficient supporting cells—it is not automatically a negative result.

`INFERENCE_target_guide_pairing_strategy` selects the local candidate pairs; `INFERENCE_max_target_distance_bp` defaults to 1 Mb. `QC_min_cells_per_gene` is a **fraction** of retained cells (default 0.05), not a count, and acts after cell QC. `QC_target_gene_rescue` can retain sufficiently observed intended target genes under explicit support thresholds. Changes to these settings change the tested universe and FDR family; record them before comparing methods or published results.

The shared covariate design includes log RNA depth and, by default, batch, detected genes and log guide UMIs. Optional mitochondrial fraction, ribosomal fraction, and assigned guide count are configurable. `cov_*` samplesheet columns are measurement-set covariates; `external_*` continuation columns are per-cell numeric covariates. Both methods receive validated requested terms. A set-level `cov_day` is typically collinear with the full measurement-set batch effect, so the pipeline rejects that design rather than silently dropping it. See [covariates](samplesheet_covariates.md) and [all parameters](parameters.md).

`INFERENCE_control_group='auto'` resolves from the declared MOI: low MOI uses non-targeting control cells; high MOI uses the complement of each perturbation. Method-specific overrides exist for sensitivity analyses and are logged when they make SCEPTRE and PerTurbo contrasts differ. Control guide annotation therefore matters biologically, not just cosmetically.

The shared resolver also checks the observed assignments: a singleton design
with control-only cells can resolve to low MOI even when high MOI was declared.
SCEPTRE's control-only cell count is computed from the input assignment matrix
and guide IDs, independently of its package-internal matrix representation.
Cells carrying targeting guides alongside controls, and cells with no guide,
are excluded from that count. This supports SCEPTRE versions that store the
internal guide matrix as sparse blocks, without changing the requested
contrast, assignments or tested pairs.

When `ENABLE_NTC_CALIBRATION=true` and both methods are run, a separate non-targeting pseudo-element pass evaluates empirical behavior without changing the main result tables or BH families. Its QQ, p-value, method-agreement and direct-target-versus-control precision–recall figures are published under `pipeline_outputs/calibration/` and shown in the live dashboard's Inference block. A skipped calibration reports why; it must not be interpreted as evidence of calibrated p-values. Compare methods on the same tested pair set and report the threshold and universe alongside AUROC/AUPRC or hit counts.


## Balanced raw-count reruns

The optional `BALANCED_source_mudata` route supports explicit cell/RNA-count
subsampling without weakening checkpoint continuation's raw-count immutability
checks. Set `ANALYSIS_mode=full`, default inference/pairing, and provide
`BALANCED_parent_inference`, `BALANCED_dataset`, positive `BALANCED_cells` and
`BALANCED_rna_umis`, and `BALANCED_seed`. The source contains the full-transcriptome
raw RNA and assigned guides; the parent inference object defines eligible cells
and the frozen analysis gene universe. Guide counts and calls must match.

Cells are selected deterministically without replacement; RNA is thinned with a
seeded joint multivariate hypergeometric draw before restricting genes. Guides
are unchanged. Both inference backends, calibration, evaluation, visualization-only
normalization, TF benchmarking and dashboards are recomputed. No old p-values
are reused. Targets must be attainable from retained cells and raw counts; a
stricter initial UMI floor may require a smaller common target. The UMI floor is
`QC_min_counts_per_cell` and applies before this optional thinning. It is not
reapplied to the downsampled counts. Inherited HTO/doublet eligibility is recorded;
subsampling does not recall doublets or guide assignment.

A single-chunk SCEPTRE handoff is copied through a temporary file, SHA-256 checked,
and atomically published. Optional additional-QC inputs have distinct staging
names, so inherited/skipped stages can safely use the same empty placeholder.
These compatibility changes preserve raw counts, assignments and tested pairs.
