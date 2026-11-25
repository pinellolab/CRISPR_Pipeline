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

## Experimental: high-MOI assignment with single-guide inference

`INFERENCE_single_guide_cells_only` defaults to `false`. Set it to `true` to
run low-MOI inference on the cells of a high-MOI screen that carry one guide,
selected at the inference boundary from the existing binary `guide_assignment`
layer and the raw guide counts, over every guide and before local/cis feature
subsetting. A cell is kept when it has **exactly one positive assignment
call**, that guide carries at least `INFERENCE_single_guide_min_umis` raw UMIs
(default 5) and at least `INFERENCE_single_guide_min_ratio` times the UMIs of
the strongest other guide in the cell (default 3). Raw guide UMIs are never
used as assignments; zero-guide and multi-guide cells, cells under the UMI
floor and cells with a competing guide are removed only here. RNA counts,
guide counts, guide identities and existing assignment calls are preserved;
no reassignment, guide collapsing, RNA thinning or new upstream QC occurs. The
selected cohort is passed to **both** SCEPTRE and PerTurbo and to calibration.

The two thresholds exist because the mixture model and a share rule fail in
opposite places. The mixture model scores each guide against its own
background, so it rejects a real second infection even when that guide is a
small fraction of a deep top guide, but on an almost empty guide profile a
single call can be noise. A share of total guide UMIs cannot tell a second
real guide from the diffuse one-UMI background of a low-depth guide library.
The ratio to the next-strongest guide is depth-independent: fourteen one-UMI
background guides under a 22-UMI top guide pass it, a 150-UMI second guide
under a 600-UMI top guide does not. Setting both thresholds to 0 reduces the
selection to "exactly one assigned guide".

Keep `Multiplicity_of_infection = 'high'` and `GUIDE_ASSIGNMENT_method = 'sceptre'`
for upstream mixture assignment. The selected inference cohort is treated as
low MOI independently of that upstream parameter. Set
`INFERENCE_control_group = 'nt_cells'` to compare with non-targeting singleton
cells for both methods; this fails if there are no control cells. Avoid stale
method-specific overrides. The default `auto` resolves from the selected cohort;
`complement` remains an explicit alternative contrast.

The pipeline publishes `inference_cell_filter/inference_single_guide_filter.json`
(retained cells, cells removed at each of the four conditions, the thresholds,
the median assigned-guide UMIs and ratio of the kept cells, and signatures of the
preserved guide matrices) and `inference_single_guide_cells.tsv`. Selection fails
if assignments are absent, nonbinary, misaligned, or no cell passes. The default
path does not add a filter process or change existing results. Enabling the
option requires fresh pair preparation, inference, calibration, evaluation and
dashboards; preprocessing and high-MOI assignments can be reused when their
exact final cell eligibility and raw counts match. Do not reuse an earlier
concatenation that precedes the parent's final QC.

A single-guide cohort can have zero control-only cells. If the resolved CRT
pool is all-cells and no control guides are assigned, PerTurbo fits its
baseline on all analysed cells without loading a separate control matrix.

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

IGVF packaging captures the requested output mode at workflow initialization
and checks the workflow status in the module binding at completion.

When the configured cell filter retains genuine zero-guide cells, both inference
backends use `log(max(total_guide_umis, 1))` for guide-depth adjustment. Negative
or non-finite guide counts remain errors; RNA depth and detected-gene counts
must be positive. This does not assign a guide to a zero-guide cell or relax
perturbation eligibility. Raw counts and assignment layers remain unchanged.
