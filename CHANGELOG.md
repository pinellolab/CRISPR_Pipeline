# nf-core/crispr: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - unreleased

### Inference

- One covariate set for both methods, derived in one place. SCEPTRE and PerTurbo now condition on the per-cell log guide-UMI depth (`total_guide_umis`), log library size (`total_gene_umis`), log detected genes (`num_expressed_genes`) and the sequencing batch (`batch`, the samplesheet's `measurement_sets`). `bin/mudata_concat.py` derives them right after the pipeline's cell filter and before its gene filter, so the cis subset, SCEPTRE's gene chunks and both inference steps inherit one set of values; `bin/prepare_inference.py` and the PerTurbo adapter only re-derive them for a user-supplied MuData. SCEPTRE receives an explicit formula (`~ log(total_guide_umis) + log(total_gene_umis) + log(num_expressed_genes) + batch`) instead of `auto_construct_formula_object`, which dropped continuous covariates with 15 or more distinct values and, on chunked runs, measured `response_n_umis` per chunk. PerTurbo receives the same quantities as precomputed `log_` columns (it then z-scores them). `percent_mito` is no longer conditioned on: it is 0 in 97% of TAP-seq cells and its distribution is panel-dependent. Validated on the TAP-seq chr8 panel: in the head run SCEPTRE reports `SCEPTRE formula: ~ log(total_guide_umis) + log(total_gene_umis) + log(num_expressed_genes) + batch` and PerTurbo reports `Covariates: continuous ['log_total_guide_umis', 'log_total_gene_umis', 'log_num_expressed_genes']; batch batch.`, so the two methods condition on the same four terms. Agreement on cis calls (3,951 pairs) rises monotonically with the change: Jaccard 0.70 with the old covariates and SCEPTRE without batch, 0.87 with the shared set, 0.91 at the branch head, and the within-run PerTurbo/SCEPTRE log10 p correlation over those pairs goes 0.956 -> 0.979 -> 0.982. All 21 pairs the TAP-seq reference reports as significant and same-direction are recovered by both methods, same direction, in all four runs; the 5 reference pairs that appear in the tables and are not recovered are the same 5 in every run and are non-significant for SCEPTRE too.
- PerTurbo's per-pair CRT diagnostics now ride along in the result tables, under a `perturbo_crt_` prefix: `crt_low_information`, `crt_observed_nonzero`, `crt_expected_nonzero`, `crt_saddlepoint_valid`, and, from the PerTurbo release that adds the Chernoff tail fallback, `crt_tail_failure_reason`, `crt_used_chernoff`, `crt_used_conservative_one` and `crt_root_residual_null_sd`. Whichever of them the PerTurbo version emits are carried; none filters a row or enters the catalog. They exist so a call can be read against how much data it rests on: on the TAP-seq chr8 panel with the shared covariates under rc9, 834 of the 899 guide-level calls that appeared were flagged low-information and 831 had no cell in which the gene was detected, sitting at the saddlepoint's 1e-12 floor. At the branch head (rc10 plus the not-yet-released Chernoff tail guard) the global per-guide table carries the diagnostics on all 277,644 PerTurbo-tested rows of 280,160: 82,825 rows are flagged `low_information`, 7,372 took the Chernoff fallback and 7 the conservative-one bound, and none of the 69 calls at q<0.05 is flagged -- on this panel the flag would gate no call that was made. Same picture in the global per-element table (15,331 of 70,788 flagged, 0 of 45 calls) and the local per-guide table (4,250 of 15,529, 0 of 88 calls); in that run the `perturbo_crt_*` columns appeared in three of the four result tables; the per-element method merge reordered its columns through a fixed list and dropped them from `local_analysis_per_element_output.parquet`, fixed in `3cd6a78` so all four tables carry them. The four tail-guard columns (`crt_tail_failure_reason`, `crt_used_chernoff`, `crt_used_conservative_one`, `crt_root_residual_null_sd`) are carried whenever the PerTurbo version emits them, and the pipeline's pinned PerTurbo, `v2.0.0rc10`, does not yet emit them, so a stock run today carries the four low-information and saddlepoint columns only.
- Fixed: the guide-UMI covariate could be dropped from an all-cells run. The adapter's identifiability guard measured `log_total_guide_umis` over the control cells alone, whatever the CRT pool, so under `all-cells` -- the shipped high-MOI default, where PerTurbo fits its stage-one baseline on every cell -- a screen with no non-targeting guides, or too few control cells to show variance, lost the covariate while SCEPTRE kept it. The guard now measures variance over the cells the resolved pool actually fits the baseline on (`_covariate_has_baseline_variance` in `bin/perturbo_v2_pipeline_adapter.py`), and still drops the covariate under `control-anchored` or `auto` when the control cells carry no variance.
- Fixed: SCEPTRE never received the batch term. `make_model_matrix_data` in `bin/inference_sceptre.R` added an empty NA level to every factor, so the rank check in the collinearity filter failed for any complete factor and dropped it as collinear, silently, on every run to date. On the TAP-seq chr8 panel SCEPTRE fit with no batch term while PerTurbo had `--batch-covariate batch`. Any earlier SCEPTRE-vs-PerTurbo comparison on this pipeline had unequal covariate sets on that axis. Dropped covariates are now printed as they are dropped. Restoring the term costs run time: on the TAP-seq chr8 panel SCEPTRE took 1,840 s with the four shared covariates against 632 s before (about 2.9x, the batch dummies being most of it), while PerTurbo stayed at 272-298 s across the runs. It also moves SCEPTRE's own calls: cis calls at q<0.05 go from 63 without the batch term to 50-51 with it, and 48 of the head run's 51 are shared with PerTurbo's 50.
- Fixed: the detected-gene count had two names and one of them lied. `bin/create_mdata.py` renamed scanpy's `n_genes_by_counts` (genes with at least one UMI) to `n_counts`, which reads as a UMI total, and renamed `n_genes` to `num_expressed_genes`, but `n_genes` only exists when `sc.pp.filter_cells(min_genes=...)` runs, i.e. without the barcode-rank filter. The count is now `num_expressed_genes` on every run and `n_counts` is no longer written. `total_gene_umis` (scanpy `total_counts`) and `num_expressed_genes` are both measured over every gene in the count matrix at QC time, before any GTF restriction; on a targeted panel that includes off-panel reads, so on TAP-seq they exceed the panel-only sums by a median of 10 UMIs and 3 genes. Both methods now name it the same way in their own fitted covariate lists -- `log(num_expressed_genes)` in SCEPTRE's formula, `log_num_expressed_genes` in PerTurbo's continuous set -- as the head TAP-seq run's logs show.

- Pin PerTurbo to `v2.0.0rc10` (`sha256:69806fe5…48dc16`). The low-MOI control-anchored CRT runs the Replogle genome-wide screen at 247-311 s per 300-target chunk against 1,013-1,274 s under rc9, restoring the pre-port pace; the screen-wide propensity basis is built on the device in row blocks (about 39 s where rc9 spent about 20 min of one host core). Against the rc9 image on Replogle essential (19,463,699 pairs) no pair is lost or gained, log10 p agrees to Pearson 1.000000000, q<0.05 calls move 901,737 -> 901,739 (Jaccard 0.999998), on-target recovery and null false discoveries are identical. On the TAP-seq chr8 panel under the all-cells pool, element-level calls and the non-targeting rate are unchanged, but **guide-level** calls at q<0.05 fall from 359 to 137: the 227 rc9-only calls were one pair each on guides with a median of 9 cells, 206 of them with no cell in which the gene is detected, sitting at the saddlepoint's 1e-12 floor - a known clipped-tail artifact, not signal, whose durable fix (a strict root-residual guard with a bounded fallback) is being ported separately. Until it lands, treat guide-level TAP-seq counts as unstable in that class; `crt_low_information` marks 91% of those pairs and is now carried into the pipeline tables as `perturbo_crt_low_information`.

The pipeline's second major version, `dev` at `f60be97`, 211 non-merge commits ahead of `refs/remotes/origin/main` (v1.0.0).

PerTurbo inference no longer rests on a posterior probability. Like SCEPTRE, PerTurbo now runs a GPU-accelerated conditional randomization test (CRT) with a saddlepoint tail beside its Bayesian effect estimates, and the CRT's p-value is what the pipeline publishes. The two inference methods were brought closer together: one control-group setting, one covariate list, one Benjamini-Hochberg family for the like-for-like comparison, and one upstream cell filter so they fit on the same cells (see below). The `cis`/`trans` vocabulary became `local_analysis`/`global_analysis` in result files, `.uns` keys and `additional_qc/` paths. Two per-pair catalog tables were added, and a series of performance fixes took the downstream steps from hours to seconds. Result tables default to Parquet.

**Why this is a major version**

1. **A v1 configuration can fail at launch.** `QC_min_cells_per_gene` given as an absolute count now
   raises (`bin/preprocess_adata.py:248`, `bin/mudata_concat.py:19`); v1 accepted `10` and ignored it in
   preprocessing while `mudata_concat` read it as a fraction.
2. **Result files were renamed.** `cis_*` / `trans_*` became `local_analysis_*` / `global_analysis_*`, in
   the published tables, the `.uns` keys and the `additional_qc/` paths, so anything downstream that opens
   those names breaks.
3. **Two parameters were removed**: `INFERENCE_PERTURBO_BATCH_SIZE` and
   `INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK`.
4. **A process was removed**: `inference_perturbo_trans`. One PerTurbo fit now produces both tables.

Upgrading from v1.0.0 also changes results. Read the next section before comparing output to v1.

### Changes that affect results

**What data is analysed**

- **Cell barcode identity changed for every run, all three modalities, unconditionally.** `bin/anndata_concat.py:252` now calls `apply_batch_suffix` for each input (only whether it *also* records `corrected_barcode` is gated on `replace_barcodes`), taking the suffix from `parse_covariate.csv`'s `barcode_key`/`concat_batch` via `get_barcode_key`, then calls `concat_on_disk(..., index_unique=None)`. v1 passed `index_unique="_"`, anndata's positional `_0`, `_1`, … in sorted-basename file order.
- **The `barcode_key` itself is new.** `subworkflows/local/prepare_mapping_pipeline/main.nf:29-54` emits `barcode_key` = the measurement set when every required modality (plus hash when `ENABLE_DATA_HASHING`) shares it, else `sample_<idx>`. v1's covariate JSON carried only `batch`.
- **The cross-modality pairing rule changed with it.** `bin/create_mdata.py:338-346` pairs modalities by exact barcode string, so an RNA cell is now matched to a guide cell by "i-th channel-order batch per modality" rather than v1's "i-th sorted file per modality". This changes the barcode strings in every output and the cell set entering inference.
- **New hard failures**: a duplicate barcode within a batch or after concatenation, and an empty modality intersection, now abort instead of producing a silently wrong MuData.
- **Cell order changed for every run.** `bin/create_mdata.py::_barcode_intersection` returns a *sorted* barcode list where v1 used `list(set(...))`, i.e. hash order, which was not stable run to run. Order-dependent randomness downstream can shift with it.
- **`DEMO_MODE` is new and ships `false`**, so a default run analyses the whole samplesheet. It defaulted to `true` for part of this release's development.
- **With `DEMO_MODE = true`, one measurement set is analysed.** `bin/filter_demo_samplesheet.py:80-89` selects the single measurement set that carries every required modality and has the fewest scRNA FASTQ files (ties broken by name); `main.nf:111` hands `params.DEMO_MODE` to `PIPELINE_INITIALISATION`, which runs `FILTER_DEMO_SAMPLESHEET` on it (`subworkflows/local/utils_nfcore_crispr_pipeline/main.nf:66-75`).
- **`DEMO_MODE` means different things on the two entry paths.** That is its only call site, and it sits in the `else` branch of the `INFERENCE_input_mudata` test (`main.nf:97-104`), so the `INFERENCE_FROM_MUDATA` entry path is not demo-filtered at all; `main.nf:138` passes the flag on to `PIPELINE_COMPLETION`, which only prints the completion warning. With `ENABLE_DATA_HASHING = true`, a samplesheet whose smallest complete set has no hash rows fails demo selection outright.

**Inference**

Both methods (PerTurbo and SCEPTRE) now test conditional independence by resampling a perturbation's assignment. That test is the CRT of Candes, Fan, Janson and Lv (2018); its power under model-X is Katsevich and Ramdas (2022); its application to single-cell CRISPR screens is SCEPTRE, Barry et al. (2021, 2024); and the saddlepoint approximation that makes it affordable at screen scale is Niu et al. (arXiv:2407.08911). Barry et al. (2024) is also where the low-MOI non-targeting-cell contrast, score test statistic, and the effective-sample-size diagnostic come from. PerTurbo (preprint coming soon) implements a fast GPU-accelerated implementation of the tests described in this literature and is not intended to exactly replicate the behavior of any prior software. The full references can be found in `CITATIONS.md`.

- **PerTurbo runs a conditional randomization test**, and its p-value is the published `p_value`. `INFERENCE_PERTURBO_CRT` defaults to `true` (`nextflow.config:117`); `_convert_common_effect_columns` takes `p_value` from `crt_saddlepoint_p_value` when the run produced one.
- **With the CRT on there is no fallback.** A pair the fit produced no CRT p-value for carries a missing `p_value` rather than the Bayesian posterior probability (`bin/perturbo_v2_pipeline_adapter.py:329-356`); BH preserves the missingness, and the posterior probability stays in `perturbo_posterior_prob`. With `INFERENCE_PERTURBO_CRT = false` the old fallback chain (`crt_p_value`, `empirical_p_value`, then the posterior probability) is unchanged.
- **Control elements are tested and carry CRT p-values.** `--crt-test-control-elements` is `argparse.BooleanOptionalAction` with `default=True` (`bin/perturbo_v2_pipeline_adapter.py:814-822`) and `modules/local/inference_perturbo/main.nf` passes no override, so the flag is always appended. Their own cells sit in the pool they are resampled within, which makes them conservative. The control evaluation now scores real p-values rather than posterior probabilities.
- The code comment at `bin/perturbo_v2_pipeline_adapter.py:345-352`, which says the CRT does not test control elements, is stale and contradicts the shipped default.
- **`QC_require_assigned_guide` (new, default `true`) drops cells carrying no assigned guide** in `bin/mudata_concat.py`, the step that runs immediately after guide assignment, so the cell population both inference methods, every derived per-cell covariate and both CRT pools see is decided once. A cell with exactly one assigned guide is kept; the count comes from `guide.layers['guide_assignment']` binarized, or from `guide.X` with a warning when that layer is absent.
- The cell filter runs before the existing gene filter, so `QC_min_cells_per_gene` is now a fraction of the cells actually analysed.
- It replaces two things: `bin/prepare_inference.py`'s `targeted_cells = grna_subset.X.sum(axis=1) > 0`, which narrowed SCEPTRE's cells a second time and is gone (its gene and guide subsetting stays), and the PerTurbo adapter's own copy of the derived `log1p_total_guide_umis_centered` formula, now a single function in `bin/inference_covariates.py` that the shared preparation step runs and the adapter consumes.
- **The filter moves numbers on the shipped `default` path.** PerTurbo's fit input is still `mudata_concat` (`subworkflows/local/inference_pipeline/main.nf:116`), but that object no longer contains cells with no assigned guide, so the all-cells CRT pool (the shipped default) now means *all cells carrying an assigned guide*, and the stage-one control fit and its size factors change with it. On the TAP-seq chr8 screen those cells were 6,831 of 126,154. Set `QC_require_assigned_guide = false` to keep them.
- The gene and guide sets still differ between the methods by design: PerTurbo fits every gene so that one run yields the transcriptome-wide table. That asymmetry is unchanged.
- **The guide-assignment rate is reported from recorded counts**, not recomputed. With the filter on, counting cells with a guide in the final MuData returns 100% by construction, so `bin/mudata_concat.py` records what it saw in the MuData's `.uns` (`n_cells_before_assigned_guide_filter`, `n_cells_after_assigned_guide_filter`, `n_cells_with_assigned_guide`, `n_cells_without_assigned_guide`, `frac_cells_with_assigned_guide`, plus which matrix it counted from and whether it filtered), whether or not it filters.
- `bin/mapping_guide.py` reports those in `additional_qc/guide/guide_metrics.tsv`'s overall row, so `frac_cells_with_guide` keeps meaning the assignment rate and is no longer `n_cells_with_guide / n_cells` once cells have been removed. A new `assigned_guide_counts_source` column says which it used; per-batch rows are still counted from the object.
- The QC metrics JSON gains an `assigned_guide_filter` block carrying the unassigned count and fraction, and the dashboard's filtering waterfall gains a step for the filter. `mean_guides_per_cell` and the guides-per-cell histogram still describe the analysed cells, so the histogram's zero bin is empty by construction.
- **`INFERENCE_control_group` (default `auto`) states, once, which cells a perturbation is compared against**, and SCEPTRE now honours it. On v1 `INFERENCE_SCEPTRE_control_group` was threaded into the driver and then discarded: `refs/remotes/origin/main:bin/inference_sceptre.R:221-225` assigned `"complement"` whatever was asked for and only warned.
- **Every low-MOI SCEPTRE p-value and effect size differs from any run made before this change.** On a low-MOI analysis SCEPTRE now compares each perturbation against the non-targeting cells, the contrast SCEPTRE documents for low MOI and the one PerTurbo's control-anchored pool already used. High-MOI analyses are unaffected, and PerTurbo's behaviour is unchanged in both cases. For the old low-MOI SCEPTRE result set use `INFERENCE_control_group = 'complement'` (which also puts PerTurbo on the all-cells pool), or override SCEPTRE alone with `INFERENCE_SCEPTRE_control_group`.
- **Under that low-MOI contrast, SCEPTRE's control gRNAs are relabelled.** They are moved out of their `non-targeting|N` buckets into SCEPTRE's reserved `non-targeting` group so they can serve as the control population (`bin/inference_sceptre.R:230-247`), which changes SCEPTRE's `grna_target` column on those rows. v1 had no such step; it aborted if the exact label appeared (`refs/remotes/origin/main:bin/inference_sceptre.R:176-177`) and always used the complement contrast. The block runs only when the assignments are low-MOI, the column is converted with `as.character` first, and the run stops if any `NA` survives.
- **The shipped defaults do not take that path.** `Multiplicity_of_infection = 'high'` (`nextflow.config:47`) with `INFERENCE_control_group = 'auto'` resolves to SCEPTRE `complement` and PerTurbo `all-cells`. A low-MOI screen has to declare `Multiplicity_of_infection = 'low'` (or set `INFERENCE_control_group` outright) before the two changes above apply to it.
- **Asking for `nt_cells` on a high-MOI analysis now fails with an explicit message** (`bin/inference_sceptre.R:314-321`) instead of silently substituting `complement`. So does asking for it when no cell carries only non-targeting guides.
- **The assignments, not the samplesheet, decide the MOI in both drivers.** A screen whose cells each carry at most one gRNA is treated as low-MOI whatever `Multiplicity_of_infection` declares (`bin/inference_sceptre.R:201-205`). Under a declared low MOI with a few double-assigned cells, SCEPTRE keeps the declared setting and lets its low-MOI QC set those cells aside rather than flipping the whole object to high MOI. Override per method with `INFERENCE_PERTURBO_CRT_POOL` / `INFERENCE_SCEPTRE_control_group`, or state the contrast outright with `INFERENCE_control_group`.
- **Both methods now condition on the same covariate list**, named once in `CANONICAL_COVARIATES` in `bin/inference_covariates.py`: `log(total_guide_umis)`, `log(total_gene_umis)`, `log(num_expressed_genes)` and the `batch` column. `percent_mito` is not among them (see the Inference section above).
- **PerTurbo's covariates changed.** v1 passed `continuous_covariates_keys=["log1p_total_guide_umis_centered"]` plus `batch_key="batch"` (`refs/remotes/origin/main:bin/perturbo_inference.py:158-163`). That derived column no longer exists; PerTurbo now receives the three counts as precomputed `log_total_guide_umis`, `log_total_gene_umis` and `log_num_expressed_genes` columns on `gene.obs` (it z-scores them itself) plus `--batch-covariate batch`.
- **SCEPTRE gains the two depth terms.** Its R reader has always fed the MuData's top-level `obs` to `import_data` as `extra_covariates`, and that frame is the intersection of the gene and guide modalities' `obs` columns (`bin/create_mdata.py:377-383`, unchanged from v1), which contains `batch` but neither depth column. `bin/mudata_concat.py` now writes the agreed columns into that frame explicitly, and `bin/prepare_inference.py` rewrites them after the cis subset rebuilds the MuData.
- The library size is deliberately absent: PerTurbo takes it as an offset and SCEPTRE adds `log(response_n_umis)` itself.
- The fifteen-level limit that `sceptre:::auto_construct_formula_object` applies (`MAX_N_LEVELS_ALLOWED`) no longer bites, because the formula is written out rather than auto-constructed. `build_formula_object` in `bin/inference_sceptre.R` fits the batch term whatever its cardinality and warns when it has fifteen or more levels, so a screen with dozens of sequencing batches conditions on the same column PerTurbo does. The claim about the builder's behaviour is still a claim about the SCEPTRE package in `sjiang9/sceptre-igvf:0.2`; nothing in this repository verifies it, and nothing now depends on it.
- **The local table is Benjamini-Hochberg-corrected over the requested pairs alone**, the same family SCEPTRE corrects over, so the two methods' q-values are comparable. v1 emitted no q-values at all.
- During this release cycle that table's family also carried every non-targeting element crossed with every tested gene: on the Replogle screen 2,623,521 of 2,714,942 rows, 96.6% null by construction, against 91,421 requested pairs (`bin/perturbo_v2_pipeline_adapter.py:220-231`). Pre-release dev outputs are therefore not comparable to the shipped behaviour either. Control pairs are still fitted, tested and present in the transcriptome-wide table, which is what the control evaluation reads.
- **One PerTurbo fit now produces both the local and the global table.** `inference_perturbo_trans` is gone; the global table is every perturbation-gene pair and the local table is the requested pairs selected out of that same fit (`--pairs-to-test` selects rows of `element_effects_requested_pairs.parquet` rather than restricting the fit). On v1 the cis numbers came from a separate fit restricted to the cis genes, guides and targeted cells, so local-analysis effect sizes and q-values differ from v1's cis numbers even at identical settings.
- **Every `q_value` column is new**, and the families differ by table. v1 emitted none: `refs/remotes/origin/main:bin/merge_method_results.py` computed no q-values and `refs/remotes/origin/main:bin/inference_sceptre.R` wrote none.
- On dev, SCEPTRE's `q_value` is `p.adjust(..., "BH")` per chunk (`bin/inference_sceptre.R:40-41`) and then recomputed across the merged chunk set (`bin/merge_sceptre_chunk_results.py`); `perturbo_q_value` is corrected within each PerTurbo table by `_bh_adjust` and recomputed over the merged table in `bin/merge_method_results.py`.
- The catalogs therefore carry `perturbo_cis_q_value` (requested pairs) beside `perturbo_q_value` (transcriptome-wide). Pairing SCEPTRE's q over the requested pairs with PerTurbo's q over every pair in the screen put a far stricter correction next to a looser one under symmetric names; `perturbo_cis_q_value` is the like-for-like number.
- **PerTurbo's SVI defaults are 500 steps per stage at learning rate 0.01** (`INFERENCE_PERTURBO_NUM_STEPS_CONTROL`, `_NUM_STEPS_BETAS`, `_STEP_SIZE`). Per the in-config note, 500 steps at 0.01 match 2,500 at 0.003, so effect-size magnitudes from any earlier 300-step run are not comparable to these. Minibatching is off (`--batch-size 0`, full batch).
- **The CRT recipe is not operator-settable.** `_run_perturbo` hardcodes `--crt-mechanism propensity --crt-tail-families saddlepoint --crt-saddlepoint-only --crt-polish-baseline --crt-allow-unconverged-baseline` (`bin/perturbo_v2_pipeline_adapter.py:504-514`). The last of those makes the null-mode guard advisory rather than fatal; no parameter exposes either the tail family or that tolerance.
- **`containers.sceptre` moved from `sjiang9/sceptre-igvf:0.1` to `:0.2`**, a version bump of the engine behind every SCEPTRE p-value, q-value and effect size. Its contents are not characterised anywhere in this repository, and the tag is mutable rather than a digest.
- **`containers.base` moved from the mutable `conda-docker:latest` to `conda-docker:sha-91bc741`.** That image is the interpreter for cell calling, QC filtering, knee detection, MuData assembly and every merge, so this is a version change of the engine behind every QC threshold and every table.
- It is also the one image this repository describes: `docker-images/conda-docker/nextflow.yaml` adds `pyarrow`, removes `muon`, and newly pins `anndata>=0.11,<0.12` while leaving scanpy and numpy unpinned, so a pinned build freezes a solve that `:latest` would have kept rebuilding.
- **`containers.perturbo` moved from the tag `perturbo:sha-f3dc8ca` to a digest pin, `perturbo@sha256:8c5d5a00…`**, described in-config as v2.0 rc9. It is a different inference engine: the v2 CRT/saddlepoint path this release is built on. A test asserts the pin keeps the `@sha256:` form (`tests/test_nextflow_container_defaults.py:27-34`).
- **rc8 (`8b25071`) carried one change over rc7**: a control-anchored CRT no longer aborts when a target perturbation has no assigned cell, 149 of 4,120 guides on the TAP-seq chr8 screen, which had killed the guide-level fit. Those targets are dropped and reported, their rows kept with missing statistics, the count recorded in `crt_metadata.json`, and the all-cells path is byte-identical.
- **rc9 (`4047280`) is documented in PerTurbo's own `CHANGELOG.md` at its `v2.0.0rc9` tag (commit `2006ab3`)**, which is the authority for the image's internals; this repository only pins the digest. Four of its five changes can move a number:
  - batch-covariate levels are enumerated over the analysed cells rather than the control cells, so a batch level holding no control cells is no longer folded silently into the reference level;
  - the all-cells propensity CRT resamples each element only within the batch levels it actually occupies, and that protection now also covers a two-level batch covariate;
  - the all-cells refit no longer anchors the batch reference on a level its fit cells do not populate;
  - the nuisance reduction over batch levels became an indicator-matrix contraction at highest precision instead of a scatter-add.
- The fifth rc9 change is additive: every CRT pair now reports how much data its p-value rests on, and a run with `--batch-covariate` reports how its control cells sit across the batch levels.
- On the TAP-seq chr8 screen, holding the pool and the gene panel fixed and changing only the image, non-targeting guides went from 21.9% to 4.7% of tests at p<0.05 against a nominal 5%, and the CRT stage fell from 1,126 s + 6,812 s to 50 s + 67 s.
- **rc9's per-pair information diagnostic does not reach the published tables.** PerTurbo writes `crt_low_information`, `crt_observed_nonzero` and `crt_expected_nonzero` into its own `element_effects.parquet`, but no file in `bin/` names any of the three, so the adapter's explicit column selection drops them; a pipeline run's per-guide and per-element outputs carry no such column (confirmed on the validated TAP-seq chr8 run). The flag marks the pairs whose p-value rests on almost no detected cells; forwarding it is a small adapter change that is not yet made.

**Pairs, panels and QC thresholds**

- **A v1 config with an absolute `QC_min_cells_per_gene` now fails the run.** `bin/preprocess_adata.py:248` raises `ValueError("Gene cell-support threshold must be a fraction in [0, 1).")`, and `bin/mudata_concat.py::resolve_min_cells` raises the same. v1 accepted the value and ignored it in preprocessing (`refs/remotes/origin/main:bin/preprocess_adata.py:190` hardcoded `min_cells=10`) while `mudata_concat` used it as a fraction, so `QC_min_cells_per_gene = 10` silently deleted every gene.
- The threshold itself did not move: `detected >= floor(n*f)+1` is equivalent to v1's `detected > int(n*f)` for every fraction in `[0, 1)`.
- **`REFERENCE_restrict_genes_to_gtf` (default `false`) keeps only the genes the resolved GTF defines** (`bin/create_mdata.py::restrict_genes_to_gtf`). For a targeted assay the GTF is the amplification panel, and whole-transcriptome mapping otherwise yields thousands of near-empty off-panel genes: on the TAP-seq chr8 screen 12,623 genes for a panel of a few dozen, with 98.5% of UMIs on the panel, each tested against every guide and each widening the BH family.
- Turning it on changes the gene panel, so it changes every q-value. Point `REFERENCE_gtf_local_path` at the panel GTF; a transcriptome-wide annotation there makes the restriction a silent no-op.
- The two panel sizes quoted for that screen describe different objects: its GTF defines 72 genes, and 68 of them survived mapping and QC into the analysed MuData, which is the number the example config and README quote.
- **`TAPSEQ_QC_MODE` (default `false`) retains every observed gene** through preprocessing (`min_cells=1`) instead of the standard 10-cell prefilter, which on a small panel can delete a real target outright.
- **`QC_min_counts_per_cell`** adds a post-barcode UMI floor (default `0`, inert).
- **The tested pair set changed.** `bin/create_pairs_to_test.py` now matches guide and GTF chromosomes through `normalize_chromosome` (leading `chr` stripped, upper-cased; lines 10-17, 28, 67-69), where v1 did an exact `query('seqname == "<guide_chr>"')`. Under the `default` pairing strategy this builds every cis candidate pair, so a library whose guide coordinates and GTF disagree on the `chr` prefix goes from zero candidate pairs to a full pair list. That is the hypothesis set itself, hence every local q-value.
- **Control-element identity changed.** `bin/intended_target_key_utils.py::annotate_intended_target_groups` now preserves explicit `element_id` groups for control guides instead of always re-bucketing them as `non-targeting|N` by median guides-per-element and lexical `guide_id` order. `bin/mudata_concat.py::preserve_source_guide_metadata` (called once, from `mudata_concat` itself) restores source-only `element_id` after concatenation so it is available to later readers.
- `annotate_intended_target_groups` is called by `create_mdata`, `create_pairs_to_test`, `prepare_inference`, `collapse_guides` and the PerTurbo adapter, so on any library that ships `element_id` on control guides the composition of the control pseudo-elements changes: the rows the control evaluation scores, and the grouping a control-anchored pool is built from. Because `collapse_guides` is one of those callers, under `DUAL_GUIDE = true` it also changes which cells survive that script's `elements_per_cell <= 1` filter, i.e. the cell set entering inference.
- **Which pairs count as on-target changed.** The new `bin/inference_target_matching.py::direct_target_mask` upper-cases and strips whitespace from intended-target identifiers, strips Ensembl version suffixes (`ENSG….3` -> `ENSG…`), maps blank strings to `pd.NA` so two absent names are never read as the same gene, and matches the target against `gene_name` as well as `gene_id`. v1 used exact string equality, `results["gene_id"] == results["intended_target_name"]` (`refs/remotes/origin/main:bin/intended_target.py:279, 422`). That changes the positive class of the intended-target and control evaluations.
- **The controls evaluation changed what it scores.** `bin/evaluate_controls.py` now picks its metric columns by preference via `select_inference_columns` (`(perturbo_log2_fc, perturbo_p_value)`, then `(sceptre_log2_fc, sceptre_p_value)`, then the generic pair), drops rows without a finite p-value from *both* classes before the 1:1 prevalence match (`keep_scorable_rows`, so an untested pair cannot move the precision-recall baseline), and takes `all_targets` from `gene_id` rather than `intended_target_name`. Control AUROC/AUPRC therefore differ from v1 beyond the posterior-probability substitution described above.
- **Genomic coordinates are normalised before they become merge keys** (`bin/merge_method_results.py`). `read_csv` types a coordinate column `int64` when complete and `float64` when one value is missing, so the same position read `100` on one side and `100.0` on the other, and the outer merge of the two methods matched nothing, leaving a table in which no row carried both methods' numbers.

**Mapping and counting**

- **The guide spacer search respects the seqspec's read id.** v1's awk expression was `t[1] = t[1] || "1"`, a logical OR that evaluates to `1` for every input, so the guide search was forced onto read 1 regardless of what the seqspec said. Dev preserves the parsed read id, handles the two-field chemistry form, and exits non-zero when it cannot infer a read id rather than guessing (`modules/local/mappingGuide/main.nf:44-79`).
- This moves where the guide is searched, so it changes guide UMIs, assignments and every test on any screen whose guide feature is not on read 1. The branch activates only when `spacer_tag` is non-empty; the shipped default is `''`, so it is off unless a site config sets a tag.
- **The 10x feature-barcode shortcut is now gated on hashing.** v1 forced `CHEM=10XV3` / `WORKFLOW=kite:10xFB` whenever `is_10x3v3` was true; dev requires `is_10xv3v == "true" && enable_data_hashing != "true"` and otherwise falls through to the seqspec-derived chemistry with plain `kite` (`modules/local/mappingGuide/main.nf:32`). A run with `is_10x3v3 = true` and `ENABLE_DATA_HASHING = true` maps guides with a different chemistry and workflow than in v1, hence different guide counts. Inert at the shipped defaults (`is_10x3v3 = false`).
- **`mappingscRNA` now passes `--sum total` and an explicit `--workflow standard` to `kb count` on the default path**, where v1 passed neither (`refs/remotes/origin/main:modules/local/mappingscRNA/main.nf:28` is a bare `kb count -i … -g … --verbose -w …`), alongside the new `--workflow nac -c1/-c2` and `--mm` paths. Whether kb treats either flag as a no-op on the standard workflow is not verifiable from this repository; see open questions.
- **Run-specific configuration was committed as the shipped defaults during development**, and has been restored to neutral values before release. `input` is `null` again and `outdir` is `'./pipeline_outputs'`, so a run without `--input` fails with a missing-input error instead of silently resolving one screen's samplesheet. `spacer_tag`, `reverse_complement_guides`, `GUIDE_ASSIGNMENT_capture_method` and `QC_barcode_filter` are back at their v1 values.
- Anyone who ran `dev` between `6a65d64` and this release got the committed values: `spacer_tag = 'TAGCTCTTAAAC'` (v1: `''`; a non-empty tag deactivates the positional guide-seqspec search and scans the whole read), `reverse_complement_guides = true` (v1: `false`; `createGuideRef` then builds the reference from reverse-complemented sequences, which collapses guide counts for a library that did not need it), `QC_barcode_filter = 'knee'` (v1: `'knee2'`; the more permissive inflection, so cell calling admits more barcodes), and `GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold` / `_n_em_rep` both `'default'` instead of `0.8` / `5`.
- Those four change guide counts, hence assignments, hence every test, so a `dev` run from that window is not comparable to a release run. The SCEPTRE assignment defaults (`'default'` rather than `0.8` / `5`) are unchanged by this restoration. Set library-specific values in a site config; `nextflow_tapseq.config` and `nextflow_cc.config` show the shape.
- **`GUIDE_ASSIGNMENT_capture_method` held `'direct-capture'` during development and is back at v1's `'CROP-seq'` at release.** It does not move guide counts at the shipped defaults either way. It is consumed in three places: `bin/create_mdata.py:265` stores it in `guide.uns["capture_method"]`, `bin/qc_metrics_json.py:334` records it in `pipeline_qc_metrics.json`, and `modules/local/guide_assignment_cleanser/main.nf:19` passes it to cleanser as `--${capture_method}`.
- Only the last changes an assignment, and only under `GUIDE_ASSIGNMENT_method = 'cleanser'`. With the shipped `'sceptre'` default, `guide_assignment_cleanser` never runs (`subworkflows/local/guide_assignment_pipeline/main.nf:14`), so the flip changes no assignment and no test, but it does change one published artifact's contents.
- **`filtered_anndata.h5ad` counts are integers, not float32.** `bin/preprocess_adata.py` dropped `adata_rna.X = adata_rna.X.astype(np.float32)` in favour of `normalize_sparse_index_dtypes`, which normalises only the sparse index arrays and deliberately preserves the compact integer count dtype.
- **A base-editing guide-mapping path is available** (`is_BaseEditing`, default `false`; `BASEEDITING_method` `'legacy'` or `'flash'`), with matching and counting conventions unlike `mappingGuide`'s kallisto/kite pseudoalignment.
- The legacy mapper (`bin/base_editing_mapping.py`) uses CRISPR-Correct ambiguity-spreading at a hardcoded Hamming threshold of 6 with reverse-complemented protospacers (`revcomp_protospacer=True`, lines 164-165) and ignores `spacer_tag` and `reverse_complement_guides`; `modules/local/mappingGuideBaseEditing/main.nf` passes neither, and nothing exposes the threshold on that path.
- The FLASH mapper takes `--tolerance` (`BASEEDITING_FLASH_tolerance`, default 5), `--reverse_complement_guides` (the shipped global, currently `true`, applied to the guide spacers at `bin/flash_base_editing_mapping.py:753-763`) and `--spacer_tag`. Nothing in the repository compares either mapper with `mappingGuide` on the same input.

### `Added`

- **`INFERENCE_control_group`** (default `auto`) names the cells a perturbation is compared against, in SCEPTRE's vocabulary, and the pipeline translates it for PerTurbo: `nt_cells` -> `--crt-pool control-anchored`, `complement` -> `--crt-pool all-cells`, `auto` -> whichever `Multiplicity_of_infection` implies (`low` -> control-anchored / `nt_cells`, `high` -> all-cells / `complement`). The mapping lives in one place (`modules/local/control_group`, a Groovy function library with no process in it, mirrored by `bin/control_group.py`), resolves once in the inference subworkflow, and logs one line naming the declared MOI, the setting, what each method resolved to, and why.
- **Control-group provenance beside the results**: `control_group_resolution.json`, a `pipeline_control_group` block in PerTurbo's `crt_metadata.json`, and `sceptre_control_group.json` from the SCEPTRE driver.
- **`catalog_per_element_output.*` and `catalog_per_guide_output.*`**, one row per (element, gene) and per (guide, gene), merging both methods' metrics with element/guide annotation and `nPerturbedCells`. Built by two new processes, `buildCatalogElement` and `buildCatalogGuide`, kept separate from `mergeMudata` so each large table caches on its own.
- **`control_batch_composition.tsv` and `control_batch_composition_note.txt`** (`854fc0c`), written by `evaluate_controls.py` into the plots directory: cells and control-only cells per `obs["batch"]` level, each level's control share and its share of all control cells.
- The note names the level the control cells concentrate in and raises a WARNING when one level holds more than 90% of them while holding less than half of the screen's cells. Both entry points write it before the result-table check, and a failure in it is caught and printed rather than taking the metrics down. Reporting only: no metric, curve or plot changes.
- **CRT memory knobs**: `INFERENCE_PERTURBO_CRT_GENE_CHUNK_SIZE` (500) and `INFERENCE_PERTURBO_CRT_MAX_GATHER_GIB` (8), separate from stage two's `INFERENCE_PERTURBO_GENE_CHUNK_SIZE` (0 = whole panel on device). Without them the test gathers every perturbation's genes at once and a transcriptome-wide screen exhausts a 40 GB card inside the saddlepoint fit. The in-config note states none of the three changes a p-value.
- **`INFERENCE_input_mudata`** reruns inference from an existing post-guide-assignment MuData, skipping mapping, QC, guide assignment and the dashboard. `-entry INFERENCE_FROM_MUDATA` invokes it explicitly, but `-entry` is not required: the default `workflow {}` itself branches on the parameter (`main.nf:96-99`), so any non-empty value redirects an ordinary `nextflow run main.nf`. Leave it at `''` for a full run. The entrypoint requires `INFERENCE_method = 'default'` and `INFERENCE_target_guide_pairing_strategy = 'default'` (`main.nf:53-61`).
- **`REFERENCE_restrict_genes_to_gtf`** and **`TAPSEQ_QC_MODE`** for targeted assays, and **`QC_min_counts_per_cell`** as a post-barcode UMI floor.
- **`nextflow_tapseq.config`**, an example site config for a targeted screen. Eight substantive deltas from `nextflow.config`: four that define the panel (`REFERENCE_restrict_genes_to_gtf`, `TAPSEQ_QC_MODE`, `QC_min_cells_per_gene`, `QC_min_genes_per_cell`) and four that describe the library chemistry (`GUIDE_ASSIGNMENT_capture_method`, `spacer_tag`, `reverse_complement_guides`, `QC_barcode_filter`), plus `input`/`outdir` reset to neutral placeholders (`null` and `./pipeline_outputs`, undoing `nextflow.config`'s committed run paths) and the `max_cpus`/`max_memory` ceilings every site has to set. `tests/test_site_config_examples.py` exempts all four of the latter as boilerplate.
- The tapseq config also carries a comment recording why `auto` was the wrong control group on that screen: all 30 non-targeting guides' cells sat in one 10x lane, 2,033 of 2,049 control-only cells.
- **`DEMO_MODE`** (default `true`) and `FILTER_DEMO_SAMPLESHEET`: trims the samplesheet to one measurement set, the complete set with the fewest scRNA FASTQ files, warns at start and completion, and publishes `DEMO_MODE_WARNING.txt` plus the filtered samplesheet. **An unedited config analyses that one set rather than the full samplesheet**; set `DEMO_MODE = false` for a real run. `PIPELINE_INITIALISATION` also passes an `enable_hashing` flag through to `bin/filter_demo_samplesheet.py` as `require_hash`, so with `DEMO_MODE = true` and `ENABLE_DATA_HASHING = true` a samplesheet whose smallest complete set has no hash rows fails demo selection outright.
- **Base-editing guide mapping**: `mappingGuideBaseEditing` (`bin/base_editing_mapping.py`) and `mappingGuideBaseEditingFlash` (`bin/flash_base_editing_mapping.py`), selected by `is_BaseEditing` and `BASEEDITING_method`, in a new `containers.bediting` image, with `BASEEDITING_FLASH_*` knobs for tolerance, guide length, device and chunk sizes.
- **CC-Perturb-seq mapping support**: `scrna_workflow` (`standard`/`nac`), `use_multimapping`, `replace_barcodes` and `bc_replacement_file`, none of which existed in v1. With `use_multimapping = true`, kb's fractional counts are rounded in `bin/anndata_concat.py` on both workflows; the `nac` path sums `mature + nascent + ambiguous` in float32 and, with `mm = false`, leaves `.X` fractional. `downloadReference` extracts `cdna.txt` and `nascent_index.txt` from the IGVF tarball for the `nac` workflow.
- **`conf/provenance.config`**, `includeConfig`'d unconditionally from line 1 of both `nextflow.config` and `nextflow_tapseq.config`. It adds a `manifest {}` block (`name = 'pinellolab/CRISPR_Pipeline'`, `author`, `homePage`, `description`, `mainScript`, `defaultBranch = 'dev'`). v1 had no manifest at all, so `workflow.manifest.*` goes from undefined to populated and the run banner and provenance output change.
- **`pipeline_info/` reproducibility files**: the last config *file* Nextflow loaded, copied verbatim to `pipeline_info/nextflow.config`, plus `nextflow.log`, `pipeline_manifest.config` (git commit, branch, remote, command line, Nextflow version, run name, success flag) and the original samplesheet. `subworkflows/nf-core/utils_nextflow_pipeline/main.nf:88-97` takes `workflow.configFiles[-1]`, so with `-c nextflow_tapseq.config` the published file is the tapseq delta, not a resolved or merged view of the run's configuration. Samplesheets may be tab- or comma-delimited and are validated up front. The timestamped `params_*.json` in that directory is v1 behaviour, unchanged.
- **`pipeline_qc_metrics.json`**, a documented machine-readable QC payload with descriptions and units, published beside `pipeline_dashboard.tar.gz`.
- **Optional, fail-open telemetry**, off by default and enabled with `-c conf/axiom.config` or `-c conf/wandb.config`. Thirteen `AXIOM_*` parameters (`nextflow.config:32-44`: enabled, dataset, token *environment variable name* (never a credential in config), ingest/API URLs, dashboard flag, two payload caps, poll and heartbeat intervals, run id/name, event file) drive `onComplete`/`onError` lifecycle events, and a W&B HTML monitor renders a pipeline dashboard. Both also switch on Nextflow's native trace (and Axiom the DAG). A failure in either is contained and does not fail the run.
- **`INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT`** (default `parquet`) selects the serialization of **every** published result table, despite its name: local-analysis, global-analysis and both catalogs (`modules/local/mergeMudata/main.nf:37` -> `bin/merge_local_global_results.py:113-115`). `inference_perturbo`'s intermediate local tables stay `.tsv.gz` regardless.
- **`INFERENCE_PERTURBO_WRITE_MUDATA`** and **`INFERENCE_SCEPTRE_WRITE_MUDATA`** (both `false`) turn the intermediate MuData writes back on; **`INFERENCE_SCEPTRE_FORCE_CHUNK`** chunks SCEPTRE even when the auto heuristic says it is unnecessary.
- **Remaining PerTurbo v2 knobs**: `INFERENCE_PERTURBO_DEVICE`, `_MAX_CHUNK_CELLS`, `_LOCAL_PARALLEL_FITS`, `_ELEMENT_GPU`, `_GUIDE_GPU`, `_JAX_CACHE_DIR`, `_SAVE_MODEL_PARAMS`, `_SIZE_FACTOR_MODE`, `_LIKELIHOOD`, `_PRIOR`.
- **SCEPTRE `fold_change` and `se_fold_change`** columns at both guide and element level (Tim Barry), and `sceptre_negLog10p` / `perturbo_negLog10p` in the analysis tables.
- **`seqSpecCheck` now publishes** to `<outdir>/pipeline_outputs/seqspeccheck`; in v1 the process had `debug true` and no `publishDir`.
- **`tf_benchmark` gained a `demo_mode` input.** `workflows/crispr_pipeline/main.nf` passes `params.DEMO_MODE` at both call sites, which appends `--demo-mode`; in that mode `bin/tf_benchmark.py` writes `DEMO_MODE_WARNING.txt` into `benchmark_output/` and overprints the figure "DEMO MODE — PRE-RUN ONLY (NOT FINAL RESULTS)" in dark red.
- **`bin/tf_benchmark_utils.py`** (new) aliases `perturbo_p_value`/`perturbo_log2_fc` onto the canonical names and raises if no p-value column survives; the benchmark had to be adapted to the column renaming.
- **Published docs (new files)**: `docs/index.html`, `docs/parameters.html`, `docs/parameters_nfcore_style.html`, `docs/mudata_schema.md`, `docs/mudata_schema_catalog.tsv`, `docs/qc_outputs.md`, `docs/qc_outputs_flat.json`, `docs/qc_outputs_flat.tsv`, `docs/filtering_workflow.{md,html}`, and `outputs/per_element_output_field_annotations.xlsx`. `docs/README.md`, `docs/usage.md` and `docs/output.md` already existed in v1 and were updated. `AGENTS.md` was added at the repo root.

### `Changed`

- **`containers.base` is pinned to `:sha-91bc741`** (`:latest` on GHCR predates the recipe change below), and new `containers.aria2` and `containers.bediting` entries were added; `containers.gmmdemux` and `containers.cleanser` are unchanged from v1.
- **`downloadReference` moved to the aria2 image.** The base image already shipped `aria2` (`docker-images/conda-docker/nextflow.yaml`, unchanged line), so the move is not about gaining the tool. `sceptre_chunk_prepare|sceptre_chunk_merge` moved from `containers.base` to `containers.sceptre`.
- **The global `process.container` default is now a closure**, `container = { params.containers.base }`, and `tests/test_nextflow_container_defaults.py::test_every_process_has_a_base_container_fallback` asserts that exact spelling stays in the global block. No process can run in the host software environment (already true in v1, whose global block set `container = params.containers.base`), and a `-c` override of `params.containers.base` now reaches it.
- **`-c` overrides of `params.containers.*` now reach the process, with one exception.** All but one of `nextflow.config`'s seventeen process-scope container assignments, and all eight in `nextflow_cc.config`, are closures rather than eagerly-evaluated bare references. Seven of v1's eight were converted; the nine added in this release were written as closures. The exception is `withName: 'demultiplex' { container = params.containers.gmmdemux }` (`nextflow.config:550`, inside the second injected `process {}` block), still bare, so a `-c` override of that one entry does not take effect.
- **`singularity.pullTimeout = '90 min'` added, and the redundant `singularity` profile removed.** Singularity was already enabled unconditionally at the top level in v1 (`refs/remotes/origin/main:nextflow.config:284-288`, the same `enabled` / `autoMounts` / `runOptions = '--nv'`), alongside a `singularity` profile that re-set `enabled` and `autoMounts` and additionally disabled conda, docker, podman, shifter, charliecloud and apptainer (`refs/remotes/origin/main:nextflow.config:181-190`); it never set `runOptions`.
- That profile is gone; the top-level block, unchanged from v1 apart from the timeout, is what applies. The ~3.6 GB PerTurbo image can exceed Nextflow's 20-minute SIF-conversion default on a network filesystem. v1's `nextflow.config` is brace-unbalanced, which changes how those scopes resolve; see open questions.
- **Result tables were renamed**, both as files and as MuData `.uns` keys: `cis_*` -> `local_analysis_*`, `trans_*` -> `global_analysis_*`. The four old `cis_*`/`trans_*` keys are now deleted from the output MuData too, alongside the generic `per_guide_results` / `per_element_results` that v1 already removed (`bin/merge_local_global_results.py:104-111`; `refs/remotes/origin/main:bin/merge_cis_trans_results.py:40-43`). `additional_qc/trans/` became `additional_qc/global_analysis/` and its plots and tables renamed with it (`bin/trans.py` gained `--prefix`, default `global_analysis`).
- The rename covers result files, `.uns` keys and `additional_qc/` paths, **not** column names: the catalogs' like-for-like PerTurbo block keeps the `perturbo_cis_*` prefix (`bin/build_catalog_per_element_output.py:31-36`). Anything reading a `cis_`/`trans_` file or key needs updating.
- **Intermediate MuData writes are off by default.** `inference_perturbo`, each `inference_sceptre` chunk and `mergedResults` no longer write `inference_mudata.h5mu`; nothing downstream read them, and `mergeMudata` assembles the published MuData from the result tables. On a Replogle-scale screen that is tens of gigabytes of writes per run removed. `INFERENCE_PERTURBO_WRITE_MUDATA` / `INFERENCE_SCEPTRE_WRITE_MUDATA` restore them.
- **Publishing is opt-in**: the default `publishDir` is `enabled: false`, so only processes declaring their own `publishDir` publish, and `pipeline_outputs/` no longer fills with full-MuData pass-through copies.
- **Library thread pools are capped inside every task** via the `env {}` scope (`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_MAX_THREADS`, `POLARS_MAX_THREADS` = 1) in all three of `nextflow.config`, `nextflow_cc.config` and `nextflow_tapseq.config`. Each example site config stands alone and must carry its own copy, which `tests/test_nextflow_container_defaults.py` enforces. Four processes opt back in inside their own script blocks: `guide_assignment_sceptre` exports `task.cpus` for BLAS, and `mergeMudata`, `buildCatalogElement` and `buildCatalogGuide` for Polars.
- Every Python process here links OpenBLAS and pyarrow and R links openblas-pthread, each sizing its pool to the whole host regardless of `task.cpus`. Measured on the TAP-seq chr8 run, `guide_assignment_sceptre` held 129 threads at `cpus = 16` with 13 tasks at once; the oversubscription cost single-threaded scripts `mapping_gene.py` 52 min inside the pipeline against 18 s alone on the same host and file, `read_h5mu` 52 min against 67 s, and the downstream QC and catalog tasks hours for minutes of work.
- **Downstream steps read result tables and specific modalities instead of whole `.uns` loads.** `mudata.read_h5mu` loads `.uns` eagerly: `backed` defers only `.X`, and `.uns` was 6.35 of 6.41 GB on the measured screen, 4.99 GB of it a 52,006,760-row guide table.
- Measured per commit on that screen: `evaluation_controls` 2.96 h / 27.3 GB -> 7.2 s / 2.0 GB (`e5c7f1c`); `evaluation_plot`'s `igv.py` 1.19 h -> 39.8 s by taking a 13,140,543-row loop off `iterrows` (`9aba601`), then 39.8 s -> 11.2 s and 15.7 GB -> 2.57 GB peak RSS by not loading the tables it never plots (`12d7477`); `intended_target.py` 34.0-34.7 s -> 7.7-8.0 s and 1.67 GB -> 0.67 GB, measured on a 2M-row subset rather than the full screen (`aa7bf9a`); and for `additional_qc_plots`, per-invocation load time 67.3 s -> 0.8-7.7 s and peak RSS 15.8 GB -> 0.3-4.9 GB (`7992398`).
- `7992398`'s own commit message retracts attributing that task's ~12.25 h wall time to the repeated loads: five loads is about 5.6 minutes and the rest was compute over 52 million rows. The 12-hour figure belongs to the thread-oversubscription fix above. Several of these tasks' wall-time before/after spans both changes, since the thread caps landed in the same week; the per-script figures are the separable ones.
  - Three of the four were checked artifact-for-artifact: `evaluation_controls` emits the same four files byte for byte, `igv.py`'s six `evaluation_output` files are md5-identical across both of its commits, and `intended_target.py`'s metrics TSV, results TSV and three PNGs are byte-identical. The `additional_qc_plots` change was verified by asserting the slim readers return the same values as `read_h5mu`, column for column across anndata's encodings, rather than by comparing artifacts.
- **The catalog builders attempt a Polars streaming plan and fall back to the pandas reference path**, with the reason printed, when the container's Polars cannot run it. No before/after timing for either catalog process exists in this repository.
- `00eb027` (15 Sep) records that the fast path had been silently falling back on *every* run from 11 Sep until that commit: `3d370c8` added `perturbo_cis_*` to both catalogs' `OUTPUT_COLUMNS` and taught the pandas path to produce them by renaming, but the streaming plan then asked for `perturbo_cis_log2_fc`, a column present in neither input, and the `ColumnNotFoundError` was swallowed into a full pandas rebuild of 52,006,760 guide-gene rows on the TAP-seq screen.
- The fast path now performs the same rename, defers a q-value still needing BH to pandas (that correction needs a global sort), and carries `negLog10p` across rather than recomputing it.
- **`INFERENCE_PERTURBO_STEP_SIZE`, `_NUM_STEPS_CONTROL` and `_NUM_STEPS_BETAS`** are 0.01, 500 and 500; minibatching is off.
- **`max_cpus` is 128** (v1: 32) **and `max_memory` is declared as the string `'256.GB'`** rather than a `MemoryUnit`, coerced with `(params.max_memory as nextflow.util.MemoryUnit)` where it is compared. Under Nextflow 26 a `MemoryUnit` assigned in a params block is stored as a plain map and every `memory` closure failed with `Cannot compare ConfigMap with MemoryUnit`, killing a full run at `prepare_covariate`.
- **Dashboard assets (`css`, `js`, `svg`) resolve against `${projectDir}`** instead of the launch directory, so a run launched elsewhere finds them. (`ENCODE_BED_DIR` was already `"${projectDir}/encode_bed_files"` in v1.) The dashboard receives the resolved `params` map as a generated JSON file rather than reading `pipeline_info/`.
- **`seqSpecCheck`'s guide matching is no longer O(reads x guides).** Guides are bucketed by length and each read scanned with one sliding window per distinct length and O(1) set lookups. This also changes which hit is reported: v1 iterated guides in list order and took the first guide that matched anywhere in the read; dev takes the leftmost position in the read at which any guide matches.
- The `seqSpecCheck` module's `debug true` was replaced by a `publishDir`, a helper block v1 had duplicated twice was removed, and `MPLCONFIGDIR`/`XDG_CACHE_HOME` now point at the task work dir because the container's `$HOME` is read-only.
- The input-tolerance behaviour (gzip detected from the `0x1f8b` magic bytes, `.fastq`/`.fq`/`.fasta`/`.fa` accepted, a directory resolved to the sequence file inside it, an unreadable input skipped with a warning) was already in v1's `bin/seqSpecCheck.py`, as were the process's `(meta, reads)` tuple input, its `reads.collate(2)` script block, and the hashing branch that calls it.
- **Counts are built sparse at a narrow dtype and never densified blindly** (`bin/count_matrix_utils.py`, applied in both the kite and base-editing paths), low-cardinality `.uns` string columns are stored as categoricals, and gzip was dropped from the intermediate MuData/AnnData writes.
- **The `local` profile retries** (`errorStrategy = { task.attempt <= 3 ? 'retry' : 'terminate' }`); v1 had that line commented out there. The `slurm` profile already had it. Default process CPUs and memory scale with `task.attempt`; the resource table changed substantially, see the appendix.
- **`example_gasperini`'s guide seqspec names its modality `crispr`**, which is what `seqspec index -m crispr` asks for. As shipped it said `guide`, the parsed chemistry came out empty, and kallisto was invoked with `-x nan`.
- **`download_pipeline/` was deleted as a duplicate of the pre-existing `download_development/`** (`e4d4400`). v1 shipped both trees side by side, and the surviving download helper is byte-identical to v1 (`git diff refs/remotes/origin/main..HEAD -- download_development` is empty).

### `Fixed`

Almost everything in this section repairs code this release introduced, not defects in v1.0.0: `observed_low_moi` appears nowhere in `refs/remotes/origin/main:bin/inference_sceptre.R`, the non-targeting factor collapse belongs to the new low-MOI control-cell path, `bin/mudata_uns_io.py`, `bin/polars_compat.py`, `bin/streaming_catalog_io.py` and both catalog builders are new files, and v1's `bin/merge_mudata.py` wrote the whole object rather than copying and patching an h5mu. The last three items are the exceptions and do change v1 behaviour.

- **Every BH helper clips p-values that overshoot `[0, 1]` and reports the count and largest excursion** instead of letting `scipy.stats.false_discovery_control` reject the array, which it does if any single element is out of range. SCEPTRE's parametric fit overshoots on three of 97,786 Replogle pairs, worst 1.0058, enough to kill that run at `buildCatalogGuide` after the fit, the merge and a 26 GB MuData write had all succeeded.
- Five copies of the helper had the same guard and all five now report before clipping (`bin/merge_method_results.py:65-79`, `bin/merge_local_global_results.py:34-46`, `bin/merge_sceptre_chunk_results.py:118-131`, and both catalog builders). The PerTurbo adapter's own `_bh_adjust` (`bin/perturbo_v2_pipeline_adapter.py:65-75`) still clips silently, with no count and no message. v1 ran no BH anywhere, so no v1 run could hit this.
- **`mergeMudata`'s `.uns` patch is durable.** `shutil.copy` does not fsync, and opening the copy with h5py in `r+` failed `EIO` on the 96-byte superblock write at offset 512 every time on the affected filesystem, reproduced in isolation. The copy is now flushed and fsynced and the open is never retried past `EIO` (`bin/mudata_uns_io.py:25-36`).
- An earlier fix retried the failing open, which was worse than the bug: the retry succeeded against a file the storage had not settled and produced a MuData whose gene matrix raised `wrong B-tree signature` on a full read. Two such files exist from two runs with different checksums, so the damage was not deterministic.
- The same commit drops `HDF5_USE_FILE_LOCKING=FALSE`: the reproduction fails with locking on, so locking was never the cause.
- **The streaming merge and enrichment run on either Polars generation**, verified on 1.44.2 and 0.20.31. `LazyFrame.collect_schema()`, `sink_parquet(engine=…)`, `collect(engine="streaming")` and `maintain_order` are Polars 1.x spellings; twenty-seven call sites across `bin/fast_result_enrichment.py`, `bin/streaming_catalog_io.py` and `bin/mudata_uns_io.py` now go through `bin/polars_compat.py` and pick the equivalent behaviour on either.
- On the older `conda-docker:latest` (Polars 0.20.31, still on the cluster) the catalog fast path cannot be taken at all: `sink_parquet` there accepts `maintain_order=True` but refuses the value with "not yet supported in standard engine". That image falls back to pandas, and `tests/test_streaming_catalog_fast_path.py` skips itself on such an image. The pinned `sha-91bc741` base image takes the fast path. The comments inside `bin/polars_compat.py` that say the pipeline's base container ships 0.20 are stale.
- **The streaming fast path survives an untested pair.** Its guard tested the q column alone, so once untested pairs carried neither p nor q every such run fell back to pandas and materialised a 211M-row enriched table. The guard now tests exactly the case the fast path cannot reproduce.
- **Every SCEPTRE task failed with `object 'observed_low_moi' not found`**, computed in the importer and read in the inference function. The inference function now takes it from the object the importer built.
- **A SCEPTRE chunk no longer dies on the new non-targeting collapse.** `grna_target` arrives as a factor, and assigning the reserved `non-targeting` label (a level it does not have) yielded `NA` with only an R warning, so the run failed later in SCEPTRE's input check with a message pointing at the data; the Replogle run lost a whole chunk to it. The column is converted with `as.character` first and the run stops if any `NA` survives (`bin/inference_sceptre.R:230-247`).
- **Control gRNAs are kept out of SCEPTRE's low-MOI discovery pairs.** Under the non-targeting-cell contrast they are the reference population, SCEPTRE refuses them as discovery targets, and their buckets are collapsed into the reserved label, so a pair naming a pseudo-element would have failed at `set_analysis_parameters`.
- **`PublishDir.setMode()` crash fixed.** `mode` cannot be a closure the way `container` can, and the default `publishDir`'s `mode` stays a bare reference for that reason.
- **`PreprocessAnnData` no longer trips Scanpy on sparse index dtypes** (`concat_on_disk` can emit a matrix whose `indices` and `indptr` differ in width, which SciPy's compiled routines reject in `eliminate_zeros`), and guide annotation dtypes and explicit dual-guide control pairs survive concatenation.
- **`igv.py`, `evaluate_controls.py` and the dashboard handle mixed chromosome metadata, missing controls, dense guide assignments, and local-only PerTurbo runs** instead of failing.
- **`debug true` was removed from every process module it was on in v1**, and the four remaining `.view()` prints, unconditional in v1 (`refs/remotes/origin/main:subworkflows/local/mapping_guide_pipeline/main.nf:40`, `mapping_hashing_pipeline/main.nf:38`, `mapping_rna_pipeline/main.nf:35`, `workflows/crispr_pipeline/main.nf:61`), now sit inside `if (params.DEBUG_VAR)` guards (`DEBUG_VAR` was already `false` by default in v1, `refs/remotes/origin/main:nextflow.config:62`). Two modules added later in this branch still declare `debug true`, `modules/local/mappingGuideBaseEditing/main.nf:3` and `modules/local/mappingGuideBaseEditingFlash/main.nf:3`, so a base-editing run still streams those per-read mappers' stdout to the console.

### `Deprecated`

- **`INFERENCE_SCEPTRE_control_group` (v1's per-method setting) and the new `INFERENCE_PERTURBO_CRT_POOL`** are per-method overrides of `INFERENCE_control_group`. Each one's shipped default (`complement`, `from-moi`) counts as unset; any other value wins for that method alone, and the pipeline logs that the two methods are then deliberately inconsistent.
- **`INFERENCE_PERTURBO_BATCH_SIZE` and `INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK` were removed**, the only two v1 parameters this release drops. Minibatching is off (full batch) and gene blocking is governed by `INFERENCE_PERTURBO_GENE_CHUNK_SIZE`, with memory otherwise bounded by `_MAX_CHUNK_CELLS`, `_CRT_GENE_CHUNK_SIZE` and `_CRT_MAX_GATHER_GIB`. Separately, `inference_perturbo` passes `--perturbation-chunk-size 0` (`modules/local/inference_perturbo/main.nf:69`) and no parameter exposes it, the same treatment as the CRT recipe flags.
- **Two parameters are documented but dead.** `INFERENCE_PERTURBO_LOCAL_MAX_CHUNK_CELLS` appears in `nextflow.config:100`, `README.md`, `nextflow_schema.json` and `docs/index.html` and in no module or script. `INFERENCE_SCEPTRE_formula_object` (`nextflow.config:138`) is referenced by no `.nf`, `.R` or `.py` file: `modules/local/inference_sceptre/main.nf` writes seven values into `args.txt` and this is not one of them. `README.md:99` says as much, and it was equally dead in v1. Setting either does nothing.
- **The `test` and `singularity` profiles are gone.** `-profile test` and `-profile singularity` now fail with `Unknown configuration profile`. `conf/test.config` is still on disk but unreachable by profile, and `conf/test_private.config` (added in this release, with `tests/csv/private_samplesheet.csv`) is unreachable the same way, since no `test_private` profile was added.
- **Code in the tree that nothing invokes.** `bin/perturbo_inference.py` and `bin/perturbo_inference_chunked.py`, which v1's `inference_perturbo` and `inference_perturbo_trans` invoked, are now referenced by nothing; `bin/chunk_mudata.py` goes with them (no `.nf` named it in v1 either, but `perturbo_inference_chunked.py` imported it, so v1 ran it). `bin/merge_cis_trans_results.py` was deleted in favour of `bin/merge_local_global_results.py`.
- `modules/local/publishFiles` was already imported-and-never-invoked in v1 (`refs/remotes/origin/main:subworkflows/local/inference_pipeline/main.nf:12`); new in this release is that its `INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT` handling is dead code from birth, and it only `cp`s and renames the extension rather than converting format.
- `bin/validate_fractional_qc_params.py`, which exists to fail a run fast on a retired absolute `QC_min_cells_per_gene`, is referenced only from `tests/test_cis_pair_and_gene_filter_utils.py`; the hard failure that actually fires comes from `bin/preprocess_adata.py` and `bin/mudata_concat.py`.
- **`nextflow_schema.json` is not enforced.** `validate_params` is accepted by `PIPELINE_INITIALISATION` and not consumed there; defining it only silences the undefined-parameter warning. Eleven schema defaults disagree with `nextflow.config`, including `QC_min_genes_per_cell` (schema `800`, config `500`) and `spacer_tag` (schema `'GAGTACATGGGG'`, config `''`). See the appendix.

### `Dependencies`

- **The base image recipe changed** (`docker-images/conda-docker/nextflow.yaml`): `pyarrow` added, `anndata` pinned to `>=0.11,<0.12`, and **`muon` removed**.
- The muon removal is load-bearing: `bin/create_mdata.py` had to switch `from muon import MuData` to `from mudata import MuData`, and no file in `bin/` or `modules/` imports muon any more. Together with pyarrow this is why `containers.base` had to move off `:latest` to `:sha-91bc741`; `:latest` on GHCR predates the change and is only rebuilt by a manual `workflow_dispatch`. The pin also freezes scanpy and numpy, which are unpinned in the recipe and which `:latest` would have kept re-solving.
- **`containers.sceptre` `sjiang9/sceptre-igvf:0.1` -> `:0.2`**, `containers.perturbo` from the tag `perturbo:sha-f3dc8ca` to the digest `@sha256:8c5d5a00…` (v2.0 rc9), and two new images: `containers.aria2 = 'biasofpriene/aria2c'` and `containers.bediting = 'bioinfolucas/crispr_mapping:latest_0236'`.
- Only the perturbo entry is a digest; the other four are mutable tags. `base` is the one image built from this repository (`docker-images/conda-docker/`, the only recipe directory in the tree); `sceptre`, `aria2` and `bediting` are third-party tags whose contents nothing here characterises.

## Detailed history (dev since v1.0.0)

211 non-merge commits ahead of `refs/remotes/origin/main`, HEAD `f60be97`. Weeks with no dev-only commits are omitted; a great deal of Jan-Feb 2026 work landed on `main` itself and is therefore already released rather than part of this delta.

Eleven commit subjects in this delta also appear on `main`: `Preserve full guide metadata in create_mdata`, `seqspeccheck: fastq not zipped are again compatible`, `Filter measurement_sets to current analysis inputs`, `Fix Nextflow defaults for QC batch column and ENCODE BED path`, `Fix PerTurbo Nextflow worker interpolation`, `Run PerTurbo chunks sequentially`, `Simplify PerTurbo efficiency mode handling`, `re-enable perturbo chunking and expose CLI options`, `config explain kmnee`, `final fix for hash seqspeck check and samples concactenation while using hash` and `update documentation`. These are rebase/cherry-pick twins: the behaviour is identical on both sides, so none of them is a v1 -> v2 change. Weeks that mention them say so.

Author dates are used for the ISO weeks. They diverge badly from commit dates across the March 2026 commits, so those weeks do not always say when a change actually reached this branch.

### 2026-W38 (Sep 14 - Sep 20)

The release week: one control-group setting for both methods, GTF-restricted genes for targeted screens, a TAP-seq example site config, and the perturbo pin moving rc7 -> rc8 -> rc9. The last of the downstream performance work landed (library thread caps in the `env {}` scope, four QC/evaluation scripts no longer loading result tables they never read), and the controls evaluation began reporting where the screen's control cells sit across its batch levels.

- `854fc0c` Report where the control cells sit across the screen's batches (Logan Blaine)
- `4047280` Pin perturbo v2.0 rc9 and keep the site-config checks to our own example (Logan Blaine)
- `08f7379` Add a TAP-seq example site config (Logan Blaine)
- `7a5a6c1` One control-group setting for both SCEPTRE and PerTurbo (Logan Blaine)
- `e5c7f1c` Stop the controls evaluation reading the table it scores (Logan Blaine)
- `aa7bf9a` Normalize intended-target identifiers per value, not per row (Logan Blaine)
- `9aba601` Take igv.py's 13.1M-row loop off the row-by-row path (Logan Blaine)
- `9a588e8` Cap the math libraries' thread pools inside every task (Logan Blaine)
- `972ac52` Let a targeted screen restrict its genes to the GTF's panel (Logan Blaine)
- `8b25071` Pin perturbo v2.0 rc8 (Logan Blaine)
- `7992398` Stop the QC scripts loading result tables they do not read (Logan Blaine)
- `6c45790` Carry the library thread caps into the CC-Perturb-seq config too (Logan Blaine)
- `34995fd` Drop an unused import from the igv tests (Logan Blaine)
- `12d7477` Stop igv.py reading the result tables it never plots (Logan Blaine)
- `00eb027` Let the catalog fast path build the columns it asks for (Logan Blaine)

`4047280` is the only container pin in the delta whose commit message has no body; nothing in the tree describes what rc9 changed. `8b25071`'s body documents rc8 in full. `7992398`'s body retracts attributing the `additional_qc_plots` task's 12.25-hour wall time to the repeated `.uns` loads; `9a588e8`'s body attributes the hours-for-minutes pattern to thread oversubscription instead. `9aba601` and `12d7477` are consecutive halves of one speedup: 1.19 h -> 39.8 s, then 39.8 s -> 11.2 s. `00eb027` records that both catalogs' streaming fast path had been silently falling back to pandas on every run since 11 Sep.

### 2026-W37 (Sep 7 - Sep 13)

The busiest week of the branch, and the one that settled the statistics: the same covariate list, control-group contrast and BH family for both methods, the CRT's p-value as the published `p_value` with no posterior-probability fallback, control elements tested, and the cis-only PerTurbo process removed so one fit yields both tables. The perturbo image moved rc2 -> rc7 across the week, gigabyte-scale intermediate MuData writes were removed from three processes, the streaming merge was made to work on either Polars generation, an `EIO`-on-reopen storage bug that had silently corrupted two output MuDatas was fixed, and a fail-open Axiom/W&B telemetry stack landed.

- `ed732c6` clean up comments in nextflow.config (Logan Blaine)
- `e483194` Pin PerTurbo to the bordered-nuisance-algebra build (Logan Blaine)
- `bd84aed` Trim redundant closing sentence from inference attribution (Logan Blaine)
- `a7008e8` Clip out-of-range p-values before BH instead of failing the run (Logan Blaine)
- `99871af` Credit SCEPTRE and spaCRT in the inference adapter (Logan Blaine)
- `8061b88` Pin perturbo v2.0 rc7 (Logan Blaine)
- `5ed6449` Update the adapter tests to the behavior the adapter actually has (Logan Blaine)
- `2198ad4` Assert the perturbo pin by digest instead of a literal tag (Logan Blaine)
- `faced68` Retry the uns patch's file open past a transient storage fault (Logan Blaine)
- `fa7f9a4` Keep Polars optional in the compatibility shim (Logan Blaine)
- `c3ea351` Give the PerTurbo step a gene-block knob beside the cell chunk (Logan Blaine)
- `6c90a1a` Collapse the non-targeting groups on characters, not on a factor (Logan Blaine)
- `696efdd` Pass the join ordering keyword only where Polars accepts it (Logan Blaine)
- `5a04759` Gate the streaming sink's keywords on the Polars that runs them (Logan Blaine)
- `4929122` Make the uns patch's copy durable, and stop retrying past EIO (Logan Blaine)
- `2854138` Run PerTurbo from the v2.0.0rc6 image (Logan Blaine)
- `0f924ed` Pass the conditional randomization test its own memory knobs (Logan Blaine)
- `0084e04` Let the streaming catalog fall back instead of failing the run (Logan Blaine)
- `ddd294e` Make the streaming merge work on the container's Polars (Logan Blaine)
- `d6c0b8a` Keep the control population out of SCEPTRE's low-MOI discovery pairs (Logan Blaine)
- `d64f497` Let SCEPTRE honour a declared low MOI when a few cells carry two guides (Logan Blaine)
- `d4814f6` Run PerTurbo from the v2.0.0rc4 image (Logan Blaine)
- `cdd20d5` Stop writing a MuData per SCEPTRE chunk (Logan Blaine)
- `caec642` Default the shared covariates to the batch alone, matching the validated runs (Logan Blaine)
- `c52dd8e` Never pass a posterior probability off as a p-value (Logan Blaine)
- `97858cb` Keep the streaming merge when a pair was never tested (Logan Blaine)
- `84cd92f` Run PerTurbo from the v2.0.0rc3 image (Logan Blaine)
- `81fed20` Accept the mitochondrial fraction under its other names (Logan Blaine)
- `798f453` Stop rewriting the MuData in the inference step, and alias rather than duplicate its tables (Logan Blaine)
- `7539f2d` Run PerTurbo from the v2.0.0rc5 image (Logan Blaine)
- `6638f1b` Name the example guide seqspec's modality the way the parser asks for it (Logan Blaine)
- `5ac120c` Normalise genomic coordinates before they become merge keys (Logan Blaine)
- `4d9cc65` Correct the local table over the requested pairs alone, as SCEPTRE does (Logan Blaine)
- `4001a4e` Give the discovery-pair filter the multiplicity it asks for (Logan Blaine)
- `3d370c8` Put both PerTurbo q-values in the catalogs and name the families (Logan Blaine)
- `2913cf6` Read backed and skip the intermediate MuData in mergedResults (Logan Blaine)
- `21c9cfa` Test the control elements, so the control evaluation has p-values to score (Logan Blaine)
- `20308a6` Let the assignments decide the pool, not the declared MOI (Logan Blaine)
- `1e53d3c` Ask for the MuData the merged-results test asserts on (Logan Blaine)
- `171fb18` Keep the mitochondrial fraction in the shared covariate set (Logan Blaine)
- `12d3cc2` Keep the memory limit a string and coerce it where it is compared (Logan Blaine)
- `0988e72` bound dashboard inference artifacts (LucasSilvaFerreira)
- `06efff1` Give both methods the same covariates, defined in one place (Logan Blaine)
- `fdbe229` Stage one at 500 steps too (Logan Blaine)
- `f8ec515` add W&B QC media smoke test (LucasSilvaFerreira)
- `f4b6bd9` apply W&B trace config in live wrapper (LucasSilvaFerreira)
- `df61e87` keep lifecycle telemetry failures contained (LucasSilvaFerreira)
- `d0a99cd` Adapter: make the CRT null-mode guard advisory, as the production Gasperini runs did (Logan Blaine)
- `c341655` fix explicit Axiom ingest URL handling (LucasSilvaFerreira)
- `ac65636` design interactive W&B pipeline dashboard (LucasSilvaFerreira)
- `a922b4f` Stage-two steps default to 300, the setting of the production Gasperini runs (Logan Blaine)
- `a60005f` stream pipeline QC to Axiom dashboards (LucasSilvaFerreira)
- `a435bfb` add strict Axiom preflight commands (LucasSilvaFerreira)
- `7b0aa6a` Allow 90 minutes to pull the PerTurbo image (Logan Blaine)
- `7248327` add fail-open Axiom pipeline telemetry (LucasSilvaFerreira)
- `680222b` Point at the rc2 image (Logan Blaine)
- `63f4a2b` One PerTurbo run yields the local and global tables; the cis-only process is gone (Logan Blaine)
- `45eedfc` consolidate W&B telemetry in HTML dashboard (LucasSilvaFerreira)
- `41fa50a` expand inference and evaluation telemetry QC (LucasSilvaFerreira)
- `1ced5c6` add fail-open live W&B HTML monitor (LucasSilvaFerreira)
- `126fd94` Stage two at learning rate 0.01 with 500 steps (Logan Blaine)
- `06a4681` use Axiom-compatible dataset name (LucasSilvaFerreira)

Two commits in this week read as contradicting each other and do not. `c52dd8e` removed the posterior-probability fallback for `p_value`, and `21c9cfa`, seventeen minutes later, turned on `--crt-test-control-elements` by default, so the control rows `c52dd8e`'s comment says are left without a p-value in fact carry CRT p-values at the shipped defaults; that comment (`bin/perturbo_v2_pipeline_adapter.py:345-352`) is stale. `41fa50a`, filed as telemetry work, is also where `bin/inference_target_matching.py` entered the tree, the module that changed how intended-target identifiers are matched.

### 2026-W36 (Aug 31 - Sep 6)

MuData file handles are closed explicitly in the PerTurbo adapter, and a generated MuData schema reference (`docs/mudata_schema.md` plus a 189-row catalog TSV) was added to the docs. `e3467c6`'s subject also appears on `main`.

- `e3467c6` update documentation (LucasSilvaFerreira)
- `33d068d` close MuData handles in PerTurbo adapter (LucasSilvaFerreira)

### 2026-W35 (Aug 24 - Aug 30)

The local/cis PerTurbo path stopped needing a patched image. A purpose-built `perturbo-v2-cis-mask` image (Dockerfile, a 437-line patch and a CI workflow) was added and deleted again four days later in favour of PerTurbo's own pair-restricted inference, which is what dev ships.

- `ed31b30` use native PerTurbo pair-restricted inference (LucasSilvaFerreira)
- `603e5b6` build combined PerTurbo v2 cis-mask image (LucasSilvaFerreira)
- `0347c29` fix evaluation for local-only PerTurbo runs (LucasSilvaFerreira)

### 2026-W34 (Aug 17 - Aug 23)

Container and dashboard hygiene: every process is now required to inherit a default container (with a test asserting it), benchmark assets resolve from the pipeline checkout rather than the launch directory, and sparse index dtypes are fixed before Scanpy QC. `181e7ed` re-added both the `containers.gmmdemux` entry and the `withName: 'demultiplex'` selector that consumes it; v1 had both (`refs/remotes/origin/main:nextflow.config:79` and `:276-277`), so it restores a rebase loss and the net `gmmdemux` delta versus v1 is nil.

- `181e7ed` Add gmmdemux container to nextflow.config (Lucas Ferreira da Silva)
- `f04c557` resolve benchmark assets from pipeline checkout (LucasSilvaFerreira)
- `7e94d31` fix dashboard guide assignment summary metrics (LucasSilvaFerreira)
- `26f8eef` Fix sparse index dtypes before Scanpy QC (LucasSilvaFerreira)
- `26489f0` enforce a default container for every process (LucasSilvaFerreira)

### 2026-W33 (Aug 10 - Aug 16)

The streaming Polars catalog path arrived: `bin/fast_result_enrichment.py` and `bin/streaming_catalog_io.py`, with `buildCatalogElement` and `buildCatalogGuide` split out of `mergeMudata` as their own processes. A post-barcode RNA UMI cell filter was added (`QC_min_counts_per_cell`), local PerTurbo pairs are masked before fitting, and the PerTurbo step size was reduced. `9227f88` only reordered the params block; `is_BaseEditing` reads `false` on both sides of it, so this is not where the base-editing default was turned off (see open questions).

- `af42e8a` Optimize Parquet catalogs and stabilize mapping inputs (LucasSilvaFerreira)
- `94f62a4` Normalize streaming catalog join key types (LucasSilvaFerreira)
- `9227f88` config change the order of the parameters (Lucas Ferreira da Silva)
- `7c9442e` Accelerate global result merge with streaming Polars (LucasSilvaFerreira)
- `c642d5c` Fix dashboard asset staging (LucasSilvaFerreira)
- `5bc1a78` reduce default number of steps (Logan Blaine)
- `cb5a1a1` Add explicit post-barcode RNA UMI cell filter (LucasSilvaFerreira)
- `6778168` Mask local PerTurbo pairs before fitting (LucasSilvaFerreira)
- `66850d2` disable perturbo v2 minibatching by default (Logan Blaine)

### 2026-W32 (Aug 3 - Aug 9)

A large consolidation week, almost all of it on Aug 7. The FLASH base-editing mapper was integrated and documented across profiles, and the January side branch's positional guide-id assignment was replaced in the merge by a spacer-keyed join that raises on duplicate or missing spacers. Parquet became the default for PerTurbo results, counts are built sparse at narrow dtypes and never densified blindly, gzip was dropped from MuData writes, `seqSpecCheck`'s O(reads x guides) guide matching was replaced by length-bucketed matching and its outputs published, `debug true` and leftover `.view()` prints were removed or gated, and the config's eager-evaluation bug was fixed so a `-c` container override reaches the process.

- `f9cd358` Fix O(reads x guides) bottleneck in seqSpecCheck guide matching (Logan Blaine)
- `f47576e` Stop exposing INFERENCE_PERTURBO_PERTURBATION_CHUNK_SIZE as a pipeline param (Logan Blaine)
- `ebdfe0c` Drop gzip compression on MuData/AnnData writes (Logan Blaine)
- `d55004f` Fix config eager-evaluation bug, dedup nextflow_cc.config, opt-in publishDir (Logan Blaine)
- `d50d164` Optimize inference result table handling (Logan Blaine)
- `d2e7c64` Apply narrow-dtype sparse construction to base_editing_mapping.py too (Logan Blaine)
- `d142611` Fix PublishDir crash: mode can't be a closure like container can (Logan Blaine)
- `a3bae5a` Integrate FLASH base-editing guide mapper (LucasSilvaFerreira)
- `a02d213` Fix stale cis/trans wording in schema and docs (Logan Blaine)
- `9d743eb` Document and configure base editing across profiles (LucasSilvaFerreira)
- `912df20` Publish seqSpecCheck outputs; gzip-compress the final merged MuData (Logan Blaine)
- `8584761` Silence startup warnings, default DEBUG_VAR off (Logan Blaine)
- `81ddce8` Preserve guide IDs across FLASH concatenation (LucasSilvaFerreira)
- `783f549` Read MuData backed and patch .uns in place instead of full read+write (Logan Blaine)
- `4e76e97` Use Parquet for PerTurbo results by default (LucasSilvaFerreira)
- `4851ce2` Gate leftover .view() debug prints behind params.DEBUG_VAR (Logan Blaine)
- `38bbd33` Inherit perturbo image and chunk-size default in nextflow_cc.config (Logan Blaine)
- `26fa916` Encode low-cardinality .uns string columns as categoricals; drop gzip on final mudata (Logan Blaine)
- `26bfd7d` Use narrow-dtype sparse matrices for count data, never densify blindly (Logan Blaine)
- `1c681ce` Pin base conda-docker image to sha-91bc741 (has pyarrow) (Logan Blaine)
- `396752a` Make FLASH base-editing mapper executable (LucasSilvaFerreira)
- `90365c5` Fix named guide index in catalog merge (LucasSilvaFerreira)
- `87ecdda` Switch PerTurbo v2 intermediate tables from TSV.gz to Parquet (Logan Blaine)
- `80c1a4a` Use Parquet-capable PerTurbo image (LucasSilvaFerreira)
- `80312de` Halve PerTurbo v2 max chunk cells for 16GB GPUs (Logan Blaine)
- `6e6403f` Load Rich pretty for scvi registry display (LucasSilvaFerreira)
- `36c9102` Remove precomputed PerTurbo recovery paths (LucasSilvaFerreira)

`e848490` in W31 had already removed `debug true` from the modules that carried it; the two base-editing mappers added in this and earlier weeks still declare it at HEAD.

### 2026-W31 (Jul 27 - Aug 2)

The PerTurbo v2 adapter (`bin/perturbo_v2_pipeline_adapter.py`) landed, replacing the v1 inference scripts, and the local/global analysis output schemas were rewritten around it: this is where `cis_*`/`trans_*` became `local_analysis_*`/`global_analysis_*` and `bin/merge_cis_trans_results.py` gave way to `bin/merge_local_global_results.py`. Runtime provenance validation, a completion manifest and the `INFERENCE_FROM_MUDATA` wiring also went in, and the perturbo container moved from a pinned sha to the `:v2-dev` moving tag and back. `3f8f425` is where `create_pairs_to_test.py` gained chromosome-key normalisation, which changes the candidate pair set.

- `c6c5f98` handle mixed chromosome metadata and missing controls (LucasSilvaFerreira)
- `acc2cdc` support dense guide assignments in dashboard plots (LucasSilvaFerreira)
- `9129594` Enable inference from existing MuData input (LucasSilvaFerreira)
- `6c7de40` Preserve guide annotation dtypes on concat (LucasSilvaFerreira)
- `6b0c470` Record pipeline source manifest on completion (LucasSilvaFerreira)
- `06866aa` Preserve explicit dual-guide control pairs (LucasSilvaFerreira)
- `e848490` Remove leftover debug true from all process modules (Logan Blaine)
- `bb726c9` reuse serialized Perturbo output during recovery (LucasSilvaFerreira)
- `bb29bdb` reuse completed cis Perturbo artifacts during recovery (LucasSilvaFerreira)
- `82d6757` Use PerTurbo image with system-PATH fix (sha-2de963e) (Logan Blaine)
- `7c907e7` Track PerTurbo :v2-dev moving tag instead of a pinned sha (Logan Blaine)
- `4d43d27` Update local and global analysis output schemas (LucasSilvaFerreira)
- `3dca376` preserve parquet through final inference merge (LucasSilvaFerreira)
- `3b4151b` add runtime provenance validation for Nextflow runs (LucasSilvaFerreira)
- `1c8287e` add parquet strategy for trans Perturbo results (LucasSilvaFerreira)
- `3f8f425` updating with public runs (LucasSilvaFerreira)
- `d1f6a36` Integrate PerTurbo v2 inference adapter (Logan Blaine)
- `3308a86` Use PerTurbo image without Docker entrypoint (Logan Blaine)
- `0ac12c4` Clean up PerTurbo v2 Nextflow parameters (Logan Blaine)

### 2026-W30 (Jul 20 - Jul 26)

`6a65d64` rewrote 349 lines of `nextflow.config`, committing a specific run's configuration as the shipped defaults: `input` and `outdir` point at that run's paths, and the file gained the injected "Auto-injected by Pipeline Configurator" blocks that override the hand-written params block above them. Most of the committed run-specific defaults in the v2.0.0 notes come from here. Container image tags were also bumped.

- `6a65d64` Capture runnable demo configuration (LucasSilvaFerreira)
- `3e89106` Update pipeline container image tags (LucasSilvaFerreira)

### 2026-W28 (Jul 6 - Jul 12)

Demo mode arrived: `FILTER_DEMO_SAMPLESHEET`, `bin/filter_demo_samplesheet.py`, the `DEMO_MODE_WARNING.txt` artifact and start/completion warnings, with `DEMO_MODE` defaulting to `true`.

- `5441b12` add demo mode (LucasSilvaFerreira)
- `c7040e9` read_me change (LucasSilvaFerreira)

### 2026-W25 (Jun 15 - Jun 21)

`bin/analysis_output_formatting.py` was added, the module that shapes the analysis tables and their annotation columns, including the `fill_null("")` that `bin/inference_target_matching.py` was later written to undo, along with resource tracing config. An `igv.py` bug that prevented plot generation was fixed. `1c26db8` also left a `.nextflow.pid` file in the tree, still tracked at HEAD.

- `a00dc49` fixing igv.py bugg preventing the generation (LucasSilvaFerreira)
- `b72b7c6` adding tracing to the utilized resourcers (LucasSilvaFerreira)
- `44a2ea8` Update README with input parameters link (Lucas Ferreira da Silva)
- `1c26db8` FIXING input file problems (LucasSilvaFerreira)

### 2026-W24 (Jun 8 - Jun 14)

PerTurbo output column names were harmonized: `perturbo_fdr_log10_p_value` (added three weeks earlier) was replaced by `perturbo_q_value` and the correction moved into each output table, so the earlier column name exists in no current output. `nextflow_schema.json` was substantially expanded and published as HTML to GitHub Pages, and the dashboard and evaluation were updated for the new field names.

- `c749355` Fixing the dashboard and evaluation (new field names compatibility) (LucasSilvaFerreira)
- `c16ebd3` schema.json added and updated readme (LucasSilvaFerreira)
- `b3069d4` hto_new_fix (LucasSilvaFerreira)
- `163c1a9` schemma visualization to git pages (LucasSilvaFerreira)
- `422ba9b` small modifications (LucasSilvaFerreira)
- `baddc1f` Harmonize PerTurbo output column names (Logan Blaine)

### 2026-W23 (Jun 1 - Jun 7)

Tim Barry's SCEPTRE `fold_change` / `se_fold_change` work from the previous week was carried into the merge and catalog builders.

- `b646cd3` added tims modification (LucasSilvaFerreira)

### 2026-W22 (May 25 - May 31)

An `INFERENCE_FROM_MUDATA` entrypoint (new `INFERENCE_input_mudata` param) reruns default inference from an existing post-guide-assignment MuData, skipping mapping, guide assignment, dashboard and QC. As wired, a non-empty value also redirects the default `workflow {}`, so `-entry` is not required to take that path. SCEPTRE outputs gained `fold_change` and `se_fold_change` at both guide and element level. `56c1e49` added an empty `tim_test_file`; `1dbd8a8` removed it.

- `82ececd` add se functionality (Tim Barry)
- `56c1e49` add test file (Tim Barry)
- `1dbd8a8` rm test file (Tim Barry)
- `ccac11d` Add MuData inference entrypoint for default workflow reruns (Logan Blaine)

### 2026-W21 (May 18 - May 24)

Two new top-level deliverables: `catalog_per_element_output` (one row per element-gene pair, merging both methods' metrics) and `pipeline_qc_metrics.json` (a documented metric payload, ~456 lines of catalog with descriptions and units). The dashboard stopped reading parameters out of `pipeline_info/` and receives the resolved `params` map as a generated JSON file.

- `c812068` qc_json_file_created (LucasSilvaFerreira)
- `0b35b54` add FDR correction to PerTurbo (Logan Blaine)
- `1be0ea2` export catalog per-element output table (Logan Blaine)

### 2026-W20 (May 11 - May 17)

Barcode-replacement mode landed with a shared cross-modality `concat_batch` key and explicit barcode suffixing, and `create_mdata` gained the fatal errors on duplicate barcodes and an empty modality intersection, plus the unconditional sorted barcode intersection that changes cell order for every run, not only under barcode replacement. The guide spacer path stopped forcing read 1 and now preserves the seqspec's read id. Runs started emitting `pipeline_info/` reproducibility files, and `download_pipeline/` was deleted as a duplicate of `download_development/`. `abc7d31` also carried three undocumented default flips (`is_10x3v3`, `spacer_tag`, the GENCODE release); all three are back to their old values in `nextflow.config` at HEAD, but `nextflow_schema.json` still advertises the flipped ones.

- `e4d4400` removing redundancy between download_development and download_pipeline. download_development is now the only one and the documentation was updated (LucasSilvaFerreira)
- `abc7d31` Export run reproducibility files to pipeline_info/ (LucasSilvaFerreira)
- `2cfa671` spacer now can work with cc-perturbseq. Default spacer is not select 1:0:0 but can select between the reads witch one shold be used to map the guide (Lucas Ferreira da Silva)
- `47a36ce` removing_ run_qc_sweep_tui.py (Lucas Ferreira da Silva)
- `7fe2d61` adding replacement (LucasSilvaFerreira)

### 2026-W17 (Apr 20 - Apr 26)

Multimapping count rounding moved out of the `nac`-only branch in `bin/anndata_concat.py`, so with `use_multimapping = true` kb's fractional counts are rounded on the `standard` workflow as well, and the unconditional float32 cast of `.X` was confined to the `nac` branch. Both are intra-branch moves: v1's `bin/anndata_concat.py` has no `mm`, `nac`, rounding or float32 code at all.

- `91c3937` last_changes_aprantly_run (LucasSilvaFerreira)

### 2026-W16 (Apr 13 - Apr 19)

The largest default change of the spring: `nextflow.config` was cleaned of the duplicated keys the March rebases had left and repointed at a non-CC screen (`reverse_complement_guides` -> true, `spacer_tag` -> `'TAGCTCTTAAAC'`, `scrna_workflow` -> `'standard'`, `replace_barcodes` -> false, `QC_min_genes_per_cell` -> 800, `QC_barcode_filter` -> `'none'`), with the CC settings moved into `nextflow_cc.config`. `da1e2ef` also made both PerTurbo scripts parse again; before it they had been Python `SyntaxError`s, so PerTurbo inference could not have started from the intervening commits.

- `da1e2ef` changes to fix the inference (should we check if the chunked has bugs?) (LucasSilvaFerreira)
- `4c9ad19` remove duplicated parameter (LucasSilvaFerreira)
- `3f13eef` adding a cc config (LucasSilvaFerreira)
- `37a05ea` adding example configs (LucasSilvaFerreira)

### 2026-W13 (Mar 23 - Mar 29)

`a025a69` repaired the rebase damage in `bin/perturbo_inference_chunked.py` and removed a duplicated `inference_perturbo` call. Two `new_config` commits moved `QC_barcode_filter` back from `'knee2'` to `'knee'`. `6b4df6a` and `5063f2d` and `063a422` share their subjects with commits on `main`: the seqSpecCheck input-tolerance behaviour they describe (gzip by magic bytes, plain `.fastq`/`.fa`, directory resolution, skip-on-unreadable) is present on both sides and is not part of this delta.

- `d55dac5` new_config (LucasSilvaFerreira)
- `5cbafed` new_config (LucasSilvaFerreira)
- `6b4df6a` seqspeccheck: fastq not zipped are again compatible (Lucas Ferreira da Silva)
- `5063f2d` final fix for hash seqspeck check and samples concactenation while using hash (LucasSilvaFerreira)
- `063a422` config explain kmnee (Lucas Ferreira da Silva)
- `a025a69` merge and fixes (LucasSilvaFerreira)

### 2026-W12 (Mar 16 - Mar 22)

PerTurbo's dataloader worker count was pinned to 0: the `NXF_TASK_CPUS` / `os.cpu_count()` autodetection was deleted and both modules pass `--num_workers 0`, so data loading is single-process regardless of `task.cpus`. On Google Batch, `machineType` for five process blocks became a function of the task's own cpus/memory. The "Run PerTurbo chunks sequentially" subject is misleading: chunks were already processed one at a time; what changed is concurrency inside each chunk. Four of these five subjects also appear on `main`.

- `d08ecbf` Implement dynamic Google Batch machine type selection based on CPU and memory requirements (LucasSilvaFerreira)
- `6f9b21d` Filter measurement_sets to current analysis inputs (Alejandro Barrera)
- `de9eea8` Run PerTurbo chunks sequentially (Logan Blaine)
- `78388e7` Run PerTurbo chunks sequentially (Logan Blaine)
- `5f4d463` Run PerTurbo chunks sequentially (Logan Blaine)

### 2026-W11 (Mar 9 - Mar 15)

Three unrelated bodies of work, each re-applied several times by rebases: PerTurbo trans-chunking CLI knobs (since superseded; the params and the chunked driver no longer exist in any module), the CC-Perturb-seq work teaching mapping the kb `nac` workflow / multimapping / barcode replacement, and `QC_batch_col`. 21 of the window's 22 commits are re-applications of 6 logical changes, and each round duplicated lines rather than replacing them: at the end of the week `inference_pipeline` invoked `inference_perturbo` twice in the same branch, `nextflow.config` declared `QC_barcode_filter` twice, and `params.outdir` was first repointed at a personal GCS scratch prefix and then deleted outright. `347a1e6` ("Preserve full guide metadata in create_mdata") is a twin of `main`'s `8cd9f66`: `refs/remotes/origin/main:bin/create_mdata.py` already merges the whole guide-metadata TSV with no column whitelist and already takes `var_names` from the merged frame's `guide_id`, and the dev-vs-main diff of that file touches neither, so this is not a v1 -> v2 behaviour change.

- `fafddbc` re-enable perturbo chunking and expose CLI options (Logan Blaine)
- `ef739a9` Fix PerTurbo Nextflow worker interpolation (Logan Blaine)
- `d3a0f27` re-enable perturbo chunking and expose CLI options (Logan Blaine)
- `ccbdd2b` Fix PerTurbo Nextflow worker interpolation (Logan Blaine)
- `c405796` Update pipeline for CC-Perturb-seq data (OlgaPushkarev)
- `aa0bcdc` re-enable perturbo chunking and expose CLI options (Logan Blaine)
- `9bee5d4` Simplify PerTurbo efficiency mode handling (Logan Blaine)
- `81e376e` Fix Nextflow defaults for QC batch column and ENCODE BED path (Alejandro Barrera)
- `766f3c3` Simplify PerTurbo efficiency mode handling (Logan Blaine)
- `716ae96` Simplify PerTurbo efficiency mode handling (Logan Blaine)
- `60a45a6` Simplify PerTurbo efficiency mode handling (Logan Blaine)
- `5f408da` Simplify PerTurbo efficiency mode handling (Logan Blaine)
- `546b8b7` Fix Nextflow defaults for QC batch column and ENCODE BED path (Alejandro Barrera)
- `4f53821` Simplify PerTurbo efficiency mode handling (Logan Blaine)
- `480d26f` re-enable perturbo chunking and expose CLI options (Logan Blaine)
- `4617bdb` re-enable perturbo chunking and expose CLI options (Logan Blaine)
- `22bb7a1` Update pipeline for CC-Perturb-seq data (OlgaPushkarev)
- `21aed61` re-enable perturbo chunking and expose CLI options (Logan Blaine)
- `1e01e55` Fix Nextflow defaults for QC batch column and ENCODE BED path (Alejandro Barrera)
- `06e07af` Fix PerTurbo Nextflow worker interpolation (Logan Blaine)
- `03f7a49` Update pipeline for CC-Perturb-seq data (OlgaPushkarev)
- `347a1e6` Preserve full guide metadata in create_mdata (Logan Blaine)

### 2026-W03 (Jan 12 - Jan 18)

A 7-line edit to the base-editing mapper changing where the guide AnnData's `var['guide_id']` comes from: the guide-design TSV rather than the observed protospacer sequence. It sat on a side branch and reached dev only in the 2026-08-07 merge, where its positional assignment was replaced by a spacer-sequence-keyed join with duplicate- and missing-spacer validation, so the mechanism this commit introduced never ran on dev.

- `eddcd52` base editing modification (LucasSilvaFerreira)

### 2025-W51 (Dec 15 - Dec 21)

The base-editing path was made runnable: a `bediting` container pin and resource block, the chemistry parser's delimiters corrected to the kb convention, the guide-metadata file staged as a `path` input and read as TSV, and the output moved to `counts_unfiltered/adata.h5ad` where the concat step looks. In-pipeline downsampling went from 400k read pairs to full depth; `0cbe72a` hardcoded a 500k cap inside the script that overrode the flag, which `da76af9` removed about six minutes later. The script's standalone argparse default is still 500000, so invoking `base_editing_mapping.py` directly downsamples unless told otherwise, and the downsample is `zcat | head -n`, the first N read pairs in file order rather than a random sample.

- `da76af9` reverting the sub sample forcing (LucasSilvaFerreira)
- `0cbe72a` base editing scripts locall implementing (LucasSilvaFerreira)

### 2025-W46 (Nov 10 - Nov 16)

An alternative guide-mapping path landed: `mappingGuideBaseEditing` driving `bin/base_editing_mapping.py`, an ambiguity-aware CRISPR-Correct mapper, in place of the kallisto/kite `mappingGuide` step. `params.is_BaseEditing` arrived defaulting to **true** with no profile, schema or README override, so an unedited config routed guide mapping onto the new path; `false` is what dev ships, but no commit in this checkout records when the default was turned off (see open questions). As committed the path could not have completed: no container selector covered the new process, its chemistry parser split on inverted delimiters versus the real `seqspec index -t kb` string, and it wrote its h5ad where the concat step does not look. All three were fixed the following month.

- `baeea4d` BASE_EDITING (LucasSilvaFerreira)

## Appendix: interface changes since v1.0.0

All values read from `nextflow.config`, `nextflow_schema.json`, `modules/local/*/main.nf` and `bin/` at dev `f60be97` against `refs/remotes/origin/main`.

### Parameters added

| Parameter | Default | What it does |
| --- | --- | --- |
| `INFERENCE_control_group` | `'auto'` | The cells a perturbation is compared against, in SCEPTRE's vocabulary (`auto`\|`nt_cells`\|`complement`), resolved once and handed to **both** methods. `auto` maps from `Multiplicity_of_infection`: low -> `nt_cells` / control-anchored, high -> `complement` / all-cells. With the shipped `Multiplicity_of_infection = 'high'` this resolves to `complement` / all-cells. |
| `INFERENCE_PERTURBO_CRT` | `true` | Run the conditional randomization test beside the Bayesian effect estimates. When on, the published `p_value` is `crt_saddlepoint_p_value` and nothing is substituted when it is missing. |
| `INFERENCE_PERTURBO_CRT_POOL` | `'from-moi'` | Per-method override for PerTurbo only (new in this release; v1 had no such parameter). `from-moi` means "not set". Any other value (`auto`\|`all-cells`\|`control-anchored`) wins for PerTurbo alone. |
| `INFERENCE_PERTURBO_CRT_GENE_CHUNK_SIZE` | `500` | Genes per chunk inside the CRT/saddlepoint fit. Memory bound only; documented in-config as not changing p-values. |
| `INFERENCE_PERTURBO_CRT_MAX_GATHER_GIB` | `8` | Cap on the gather buffer the CRT allocates per chunk. Memory bound only. |
| `INFERENCE_PERTURBO_DEVICE` | `'gpu'` | `--device` to the v2 adapter. |
| `INFERENCE_PERTURBO_NUM_STEPS_CONTROL` | `500` | SVI steps for stage one (control/baseline fit). |
| `INFERENCE_PERTURBO_NUM_STEPS_BETAS` | `500` | SVI steps for stage two (effect coefficients). |
| `INFERENCE_PERTURBO_STEP_SIZE` | `0.01` | Adam learning rate for both SVI stages. In-config note: 500 beta steps at 0.01 match 2,500 at 0.003. |
| `INFERENCE_PERTURBO_MAX_CHUNK_CELLS` | `20000` | Cells per device chunk during the fit. |
| `INFERENCE_PERTURBO_GENE_CHUNK_SIZE` | `0` | Stage-two genes fitted at a time; `0` keeps the whole panel on device. |
| `INFERENCE_PERTURBO_LOCAL_MAX_CHUNK_CELLS` | `20000` | **Dead.** In `nextflow.config`, `README.md`, the schema and `docs/index.html`; referenced by no module or script. |
| `INFERENCE_PERTURBO_LOCAL_PARALLEL_FITS` | `false` | Fit element- and guide-level models concurrently on two GPUs. |
| `INFERENCE_PERTURBO_ELEMENT_GPU` | `'0'` | Device index for the element-level fit. |
| `INFERENCE_PERTURBO_GUIDE_GPU` | `'0'` | Device index for the guide-level fit. Schema says `'1'`. |
| `INFERENCE_PERTURBO_JAX_CACHE_DIR` | `'.perturbo_jax_cache'` | XLA compilation cache directory. |
| `INFERENCE_PERTURBO_SAVE_MODEL_PARAMS` | `false` | Dump fitted variational parameters into the v2 artifact dir. |
| `INFERENCE_PERTURBO_SIZE_FACTOR_MODE` | `'observed'` | Size-factor treatment in the PerTurbo likelihood. |
| `INFERENCE_PERTURBO_LIKELIHOOD` | `'negbin'` | Observation model. |
| `INFERENCE_PERTURBO_PRIOR` | `'normal'` | Prior on the effect coefficients. |
| `INFERENCE_PERTURBO_WRITE_MUDATA` | `false` | Whether `inference_perturbo` emits the intermediate `inference_mudata.h5mu`. Makes that process output optional. |
| `INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT` | `'parquet'` | `'parquet'` or `'tsv.gz'` for **every** published result table: local-analysis, global-analysis and both catalogs. Read by `inference_perturbo`, `mergeMudata`, `buildCatalogElement`, `buildCatalogGuide` (and by `publishFiles`, which is imported but never invoked). `inference_perturbo`'s intermediate local tables are hardcoded `.tsv.gz` regardless. |
| `INFERENCE_SCEPTRE_WRITE_MUDATA` | `false` | Exported as `SCEPTRE_WRITE_MUDATA`; whether each `inference_sceptre` chunk writes its own MuData. Makes that output optional. |
| `INFERENCE_SCEPTRE_FORCE_CHUNK` | `false` | `--force-chunk` to `sceptre_chunk_prepare`: chunk even when the auto heuristic says it is unnecessary. |
| `INFERENCE_input_mudata` | `''` | Path to an existing `.h5mu`. A non-empty value redirects the default `workflow {}` to `INFERENCE_FROM_MUDATA` (`main.nf:96-99`); no `-entry` needed. Requires `INFERENCE_method='default'` and `INFERENCE_target_guide_pairing_strategy='default'`. Schema says `null`. |
| `REFERENCE_restrict_genes_to_gtf` | `false` | Keep only genes the resolved GTF defines, so a targeted panel is not tested against thousands of off-panel genes. **Not present in `nextflow_schema.json`.** |
| `TAPSEQ_QC_MODE` | `false` | Retain every observed gene in `PreprocessAnnData` (`min_cells=1`) instead of the standard 10-cell prefilter. |
| `QC_min_counts_per_cell` | `0` | Post-barcode UMI floor in `PreprocessAnnData`, alongside `QC_min_genes_per_cell`. |
| `DEMO_MODE` | `false` | `FILTER_DEMO_SAMPLESHEET` trims the samplesheet to the single complete measurement set with the fewest scRNA FASTQ files, warns at start and completion, and writes `DEMO_MODE_WARNING.txt`. Also appends `--demo-mode` to `tf_benchmark`. Off by default; schema agrees. |
| `is_BaseEditing` | `false` | Route `mapping_guide_pipeline` to the base-editing mappers instead of `mappingGuide`. Opt-in as shipped; an intermediate branch state briefly defaulted it on (`baeea4d`, 2025-11-14). See open questions for when it was turned off, which no commit in this checkout records. |
| `BASEEDITING_method` | `'legacy'` | `'legacy'` -> `mappingGuideBaseEditing`, `'flash'` -> `mappingGuideBaseEditingFlash`; anything else is a hard error. |
| `BASEEDITING_FLASH_tolerance` | `5` | Mismatch tolerance for the FLASH mapper (the legacy mapper's Hamming threshold of 6 is hardcoded and unexposed). |
| `BASEEDITING_FLASH_guide_len` | `0` | Guide length; `0` infers it. |
| `BASEEDITING_FLASH_device` | `'auto'` | CPU/GPU selection. |
| `BASEEDITING_FLASH_fastq_chunk_size` | `200000` | Reads per FASTQ chunk. |
| `BASEEDITING_FLASH_gpu_read_chunk` | `4096` | Reads per GPU batch. |
| `scrna_workflow` | `'standard'` | `mappingscRNA` kb workflow; anything other than `standard` uses `nac` with `-c1 cdna` / `-c2 nascent_idx`. |
| `use_multimapping` | `false` | `--mm` to `mappingscRNA` and `anndata_concat`; the concat step then rounds kb's fractional counts on either workflow. |
| `replace_barcodes` | `false` | Substitute cell barcodes via `bc_replacement_file`; changes the counts directory to `counts_unfiltered_modified`. |
| `bc_replacement_file` | `''` | TSV of barcode replacements; staged only when `replace_barcodes` is true, otherwise `dummy_bc_replacement.txt`. |
| `validate_params` | `true` | Accepted by `PIPELINE_INITIALISATION` but not consumed there. Schema validation is not wired up; defining it only silences the undefined-parameter warning. Absent from the schema itself. |
| `AXIOM_*` (13 params) | off | `AXIOM_enabled` (`false`), `_dataset`, `_token_env` (env var *name*, not a credential), `_ingest_url`, `_api_url`, `_create_dashboard`, `_max_bytes`, `_max_event_bytes`, `_poll_seconds`, `_heartbeat_seconds`, `_run_id`, `_run_name`, `_event_file` (`nextflow.config:32-44`). Fail-open lifecycle telemetry from the `onComplete`/`onError` hooks; `conf/axiom.config` flips it on and also enables trace + DAG. |
| `WANDB_*` (7 params) | not in `nextflow.config` | `WANDB_enabled`, `_project`, `_entity`, `_token_env`, `_poll_seconds`, `_max_bytes`, `_layout`, defined only in `conf/wandb.config`, so they exist only when that file is passed with `-c`. |

### Parameters removed

The v1 -> v2 parameter removals are exactly two.

| Parameter | Old default | Replacement |
| --- | --- | --- |
| `INFERENCE_PERTURBO_BATCH_SIZE` | `4096` | None. `inference_perturbo` hardcodes `--batch-size 0` (full batch) and `--perturbation-chunk-size 0`; memory is bounded by `_MAX_CHUNK_CELLS`, `_GENE_CHUNK_SIZE`, `_CRT_GENE_CHUNK_SIZE`, `_CRT_MAX_GATHER_GIB`. |
| `INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK` | `8000` | `INFERENCE_PERTURBO_GENE_CHUNK_SIZE` (`0` = whole panel on device). The trans-only chunked driver that consumed it is called by no module. |

`--perturbation-chunk-size` is hardcoded to `0` by `modules/local/inference_perturbo/main.nf:69` and no parameter exposes it. There was never a `INFERENCE_PERTURBO_PERTURBATION_CHUNK_SIZE` in v1 to remove: it was added and unexposed again inside this cycle (`d1f6a36`, then `f47576e`).

### Parameter defaults changed

| Parameter | main | dev | Effect |
| --- | --- | --- | --- |
| `containers.perturbo` | `perturbo:sha-f3dc8ca` | `perturbo@sha256:8c5d5a00…` (in-config: v2.0 rc9) | A different inference engine: the v2 CRT/saddlepoint path. Digest pin, so it cannot drift; a test asserts the `@sha256:` form. rc8's contents are documented in `8b25071`; rc9's in PerTurbo's own CHANGELOG at its `v2.0.0rc9` tag. |
| `containers.base` | `conda-docker:latest` | `conda-docker:sha-91bc741` | Pinned so merge/catalog steps get pyarrow and `mudata` without `muon`; `:latest` on GHCR predates the recipe change. Also the interpreter for QC, cell calling and every merge, so this is an engine change, not only a packaging one. |
| `containers.sceptre` | `sjiang9/sceptre-igvf:0.1` | `sjiang9/sceptre-igvf:0.2` | New SCEPTRE image behind every SCEPTRE statistic, now also used by `sceptre_chunk_prepare`/`_merge`. Contents not characterised here. |
| `input` | `null` | `null` (unchanged at release) | Held one run's samplesheet path during development (`6a65d64`); restored before release, so a run without `--input` fails with a missing-input error as in v1. |
| `outdir` | `'./pipeline_outputs'` | `'./pipeline_outputs'` (unchanged at release) | Held one run's output prefix during development; restored before release. |
| `spacer_tag` | `''` | `''` (unchanged at release) | Held `'TAGCTCTTAAAC'` during development: a non-empty tag deactivates the positional guide-seqspec search, scans the whole read, and activates the read-id-preserving awk branch. Changes guide counts. Schema says `'GAGTACATGGGG'`. |
| `reverse_complement_guides` | `false` | `true` | `createGuideRef` builds the reference from reverse-complemented sequences. Wrong for a library that did not need it: guide counts collapse. Schema says `false`. |
| `GUIDE_ASSIGNMENT_capture_method` | `'CROP-seq'` | `'CROP-seq'` (unchanged at release) | Held `'direct-capture'` during development. Changes the assignment model **only** under `GUIDE_ASSIGNMENT_method = 'cleanser'`, which never runs at the shipped `'sceptre'` default; otherwise stored in `guide.uns["capture_method"]` and recorded in `pipeline_qc_metrics.json`. Schema says `'crop-seq'`. |
| `QC_barcode_filter` | `'knee2'` | `'knee2'` (unchanged at release) | Held `'knee'`, the more permissive inflection, in the injected block during development. Schema says `'knee2'`; v1's schema said `'stringent'` and already disagreed with v1's config. |
| `GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold` | `0.8` | `'default'` | Hands the decision to SCEPTRE's own default; changes which guide-cell assignments are called. Schema still says `0.8`. |
| `GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep` | `5` | `'default'` | Same, for EM replicates. Schema still says `5`. |
| `INFERENCE_predefined_pairs_to_test` | `null` | `'path/to/file.csv'` | Only read under the `predefined_pairs` strategy, so the `default` strategy is unaffected, but the value is no longer a usable null sentinel. Schema still says `null`. |
| `max_cpus` | `32` | `128` | Raises every `Math.min(N, params.max_cpus)` ceiling. Re-stated as `128` in both injected resource blocks. Absent from the schema. |
| `max_memory` | `256.GB` (MemoryUnit) | `'256.GB'` (String) | Same effective size; the type change is why every memory closure now casts `(params.max_memory as nextflow.util.MemoryUnit)`. Under Nextflow 26 a MemoryUnit in a params block becomes a plain map and every closure failed. The first injected block sets `'600.GB'`; the later one resets it to `'256.GB'`, and last-wins. Absent from the schema. |
| `css` / `js` / `svg` | `'assets/css'` etc. (launch-relative) | `"${projectDir}/assets/css"` etc. | Dashboard assets resolve against the pipeline checkout, so a run launched elsewhere finds them. (`ENCODE_BED_DIR` was already `"${projectDir}/encode_bed_files"` in v1 and is unchanged.) |

`is_BaseEditing` and `DEMO_MODE` are new parameters rather than changed defaults; both are listed under "Parameters added" above.

**Four injected `params {}` blocks and two injected `process {}` blocks** (`nextflow.config:475-628`) sit *after* the hand-written params block and therefore win: two headed `// --- Auto-injected by Pipeline Configurator ---` (lines 476, 559) and two headed `// --- Computational resource overrides injected by Pipeline Configurator ---` (lines 491, 573), plus `process {}` blocks at 496 and 578.

Of the parameters those blocks set, only three also have a hand-written value earlier in the file that is therefore silently discarded: `max_cpus` (line 160), `max_memory` (line 161) and `INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT` (line 114). The rest (`scrna_workflow`, `use_multimapping`, `replace_barcodes`, `bc_replacement_file`, `QC_barcode_filter`, `ENABLE_BENCHMARK`, `QC_batch_col`, `INFERENCE_input_mudata`, `ENCODE_BED_DIR`, `DEBUG_VAR`) exist *only* in the injected blocks. The overridden SCEPTRE parameters are exactly five: `INFERENCE_SCEPTRE_CHUNK_MODE`, `_GENE_CHUNK_SIZE`, `_WRITE_MUDATA` (lines 482-486) and `_MAX_MATRIX_ENTRIES`, `_FORCE_CHUNK` (lines 563-564). The hand-written `INFERENCE_SCEPTRE_control_group`, `_side`, `_grna_integration_strategy`, `_resampling_approximation`, `_resampling_mechanism` and `_formula_object` (lines 128-138) are **not** re-declared and are not discarded; editing them there works.

### Resources and machine types

| Selector | main | dev |
| --- | --- | --- |
| `process` (default) | `cpus = min(16, max_cpus)`, `memory = 50.GB * attempt`, `machineType = n2-highmem-16` on google, `container = params.containers.base` (bare), `publishDir` enabled | the default block keeps only `container = { … }` and `publishDir … enabled: false`; `cpus`, `memory` and the google `machineType` are dropped from it. The injected blocks supply `cpus = min(16 * attempt, max_cpus)` and `memory = 200.GB * attempt` (the later block `50.GB * attempt`, last-wins) but **no default `machineType`**, so on `-profile google` a process matched by no `withName` selector now gets none |
| 37 named processes (`downloadGTF\|…\|evaluation_controls`, `nextflow.config:334`) | — (no such selector; these inherited the default block's `cpus = min(16, max_cpus)`) | `container = { params.containers.base }`, `cpus = min(4 * attempt, max_cpus)`, `memory = 50.GB * attempt`, `n2-highmem-16` on google. The largest resource change in the release, cutting those processes' CPU ceiling from 16 to 4 at first attempt |
| `mappingGuide\|mappingHashing\|mappingscRNA` | `cpus = min(30, max_cpus)`, `memory = 100.GB * attempt`, `n2-highmem-32`, `disk = '200 GB'` | `cpus = 8`, `memory = 32.GB`, `n4-standard-8`, `disk = 300.GB`; injected block adds `maxRetries = 2` and `errorStrategy = { attempt <= 2 ? 'retry' : 'terminate' }`. A cut in CPU and memory, not a raise |
| `guide_assignment_cleanser` | `cpus = min(16, max_cpus)`, `memory = 100.GB * attempt`, `n2-highmem-32` on google, `scratch = true`, `disk = '500 GB'` | `cpus = min(8 * attempt, max_cpus)`, `memory = 300.GB * attempt`, `machineType = 'n2-highmem-32'` unconditionally, `scratch = '/tmp'`, `disk = '500 GB'`; injected block resets `scratch = true` |
| `guide_assignment_sceptre` | shared selector with `inference_sceptre`: `cpus = min(16, max_cpus)`, `memory = 200.GB * attempt`, `n2-highmem-64` | split out: `cpus = 16`, `memory = 128.GB`, `n2-highmem-16` on google. Exports `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS` = `task.cpus` in its own script |
| `inference_sceptre` | as above | `cpus = min(8 * attempt, max_cpus)`, `memory = 100.GB`, `n2-highmem-32` on google. This value is passed straight through as SCEPTRE's `n_processors` (`args.txt` line 7) |
| `inference_perturbo` | selector `inference_perturbo\|inference_perturbo_trans`: `cpus = min(8, max_cpus)`, `memory = 40.GB * attempt`, `machineType = 'a2-ultragpu-1g'`, `// accelerator = [ request: 1, type: 'nvidia-a100-80gb' ]` commented out | `cpus = min(8, max_cpus)`, `memory = 100.GB * attempt`, **`machineType = 'a2-highgpu-1g'`** (40 GB A100, down from 80 GB), **`accelerator = [request: 1, type: 'nvidia-tesla-a100']` active** (a different accelerator type than v1's commented-out line, not that line uncommented), `clusterOptions = '--gres=gpu:4'` on slurm. On a 40 GB card a transcriptome-wide screen needs the CRT memory knobs; `INFERENCE_PERTURBO_GENE_CHUNK_SIZE` can go back up on an 80 GB card. Neither changes a p-value: this is scheduling and host memory, i.e. whether the run completes |
| `sceptre_chunk_prepare\|sceptre_chunk_merge` | `containers.base`, `cpus = min(8 * attempt, max_cpus)`, `memory = 100.GB * attempt`, `n2-highmem-32` | `containers.sceptre`, `cpus = min(8, max_cpus)`, `memory = 100.GB` |
| `downloadReference` | `cpus = min(16, max_cpus)`, `memory = 100.GB * attempt`, `n1-highmem-32` | `containers.aria2`, `cpus = min(4 * attempt, 48)`, otherwise unchanged. (The base image already shipped `aria2`, so the image move is not about gaining the tool.) |
| `mergeMudata\|buildCatalogElement\|buildCatalogGuide` | — (no such selector) | own `process {}` block: `containers.base`, `cpus = min(16, max_cpus)`, `memory = min(200.GB, max_memory)`; each exports `POLARS_MAX_THREADS = task.cpus` in its script |
| `mappingGuideBaseEditing\|mappingGuideBaseEditingFlash` | — | `containers.bediting`, `cpus = min(8 * attempt, max_cpus)`, `memory = 200.GB * attempt`, `n2-highmem-64` on google |
| `additional_qc_plots` | inherits the default block | injected block: `cpus = min(4 * attempt, max_cpus)`, `memory = 50.GB * attempt`, `n2-highmem-16` on google |
| `env {}` scope | — | `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, `NUMEXPR_MAX_THREADS`, `POLARS_MAX_THREADS` all `'1'`, in `nextflow.config`, `nextflow_cc.config` and `nextflow_tapseq.config` |
| `local` profile | `errorStrategy` line commented out | `errorStrategy = { task.attempt <= 3 ? 'retry' : 'terminate' }` (the `slurm` profile already had it) |

### Processes added

| Process | What it does |
| --- | --- |
| `buildCatalogElement` | Runs `build_catalog_per_element_output.py` over the local + global per-element tables and the merged MuData; emits `catalog_per_element_output.{parquet\|tsv.gz}`. Wired in after `mergeMudata`, deliberately separate so each large table caches on its own. |
| `buildCatalogGuide` | The same for the per-guide catalog view; emits `catalog_per_guide_output.{parquet\|tsv.gz}`. |
| `mappingGuideBaseEditing` | Legacy base-editing guide mapper (`bin/base_editing_mapping.py`), selected by `is_BaseEditing` + `BASEEDITING_method='legacy'`; runs in `containers.bediting`. Still declares `debug true`. |
| `mappingGuideBaseEditingFlash` | FLASH base-editing guide mapper (`bin/flash_base_editing_mapping.py`), selected by `BASEEDITING_method='flash'`; also `containers.bediting`. Still declares `debug true`. |
| `filter_demo_samplesheet` / `FILTER_DEMO_SAMPLESHEET` | Runs in `PIPELINE_INITIALISATION` when `DEMO_MODE` is true; emits the reduced samplesheet plus `DEMO_MODE_WARNING.txt`, published to `pipeline_info/` and the top level of `--outdir`. Takes an `enable_hashing` flag that becomes `--require-hash`. |
| `modules/local/control_group` | A new module directory containing **no process**: a Groovy function library (`resolveControlGroupFromParams`, `resolveControlGroup`, `controlGroup*BySetting/ByMoi`) included by `inference_pipeline` and mirrored by `bin/control_group.py`. |

`tf_benchmark` gained a fourth input, `val demo_mode`, passed as `params.DEMO_MODE` at both call sites.

### Processes and profiles removed

| Removed | Note |
| --- | --- |
| `inference_perturbo_trans` | Deleted with its include. It was a second PerTurbo pass over `concat_mudata` with a dummy channel forcing it after the cis pass; `inference_perturbo` now emits both the local and the global table from one fit, run on `mudata_concat`. |
| `-profile test` | `test { includeConfig 'conf/test.config' }` is gone from the profiles block, which now holds only `local`, `slurm` and `google`. `conf/test.config` is still on disk but unreachable by profile; `-profile test` fails with `Unknown configuration profile`. `conf/test_private.config`, added in this release, has no profile either. |
| `-profile singularity` | Gone. It had re-set `singularity.enabled`/`autoMounts` and disabled conda, docker, podman, shifter, charliecloud and apptainer. Singularity was already enabled unconditionally at the top level in v1; that block now also sets `pullTimeout = '90 min'`. |

### Outputs renamed

| main | dev |
| --- | --- |
| `cis_per_element_results.tsv.gz` | `local_analysis_per_element_output.{parquet\|tsv.gz}` |
| `cis_per_guide_results.tsv.gz` | `local_analysis_per_guide_output.*` |
| `trans_per_element_results.tsv.gz` | `global_analysis_per_element_output.*` |
| `trans_per_guide_results.tsv.gz` | `global_analysis_per_guide_output.*` |
| `.uns` keys `cis_per_{guide,element}_results`, `trans_per_{guide,element}_results` | `local_analysis_per_{guide,element}_results`, `global_analysis_per_{guide,element}_results`. The old four are now deleted from the output MuData as well, alongside the generic `per_guide_results` / `per_element_results` that v1 already deleted. |
| `perturbo_cis_per_{element,guide}_output.tsv.gz` (work dir) | `perturbo_local_analysis_per_{element,guide}_output.tsv.gz`, plus new `perturbo_global_analysis_per_{element,guide}_output.*`. The `perturbo_trans_*` pair is gone with `inference_perturbo_trans`. |
| `additional_qc/trans/`, `trans_*.{tsv,png}` | `additional_qc/global_analysis/`, `global_analysis_*` (`bin/trans.py` gained `--prefix`, default `global_analysis`). |
| `trans_perturbo_{precision_recall_roc,volcano_plot,barplot_direct_vs_control}.png` | `global_analysis_perturbo_*` |

Extensions follow `INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT`, so the default local, global and catalog outputs are `.parquet`, not `.tsv.gz`. The rename does **not** extend to column names: the catalogs keep a `perturbo_cis_*` block by design.

### Outputs added

| Output | Contents |
| --- | --- |
| `catalog_per_element_output.{parquet\|tsv.gz}` | One row per (element, gene). `sceptre_{log2_fc,p_value,q_value,fc_se,negLog10p}`, `perturbo_cis_*` (the same five, over the requested-pair family; the like-for-like comparison to SCEPTRE), `perturbo_*` (the same five, transcriptome-wide), then `element_id`, `element_type`, `element_chr/start/end`, `element_name`, `guide_ids`, `num_guides`, `gene_name`, `gene_id`, `nPerturbedCells`. |
| `catalog_per_guide_output.{parquet\|tsv.gz}` | One row per (guide_id, gene_id); the same metric columns plus `guide_id`, `guide_sequence`, `guide_type`, `targeting`, `guide_chr/start/end/strand`, `pam`, `intended_target_name/_chr/_start/_end`, `gene_name`, `gene_id`, `nPerturbedCells`. |
| New columns in the local/global tables | `sceptre_negLog10p`, `perturbo_negLog10p` (floored at 1e-300); element tables gain `element_id/type/chr/start/end/name`, `guide_ids`, `num_guides`, `gene_name`, `nPerturbedCells`; guide tables gain `guide_sequence`, `guide_type`, `targeting`, `guide_chr/start/end/strand`, `pam`, `gene_name`, `nPerturbedCells`. PerTurbo rows also carry `perturbo_fc_se`, `perturbo_posterior_prob`, `perturbo_q_value`. SCEPTRE rows gain `fold_change`, `se_fold_change` and `q_value`; v1 emitted no q-values from either method. |
| `perturbo_v2_outputs/` (optional) | Raw PerTurbo v2 artifacts: `element/` and `guide/` subtrees (`element_effects.parquet`, `element_effects_requested_pairs.parquet`, `crt_metadata.json`), `element_pairs_to_test.parquet`, `guide_pairs_to_test.parquet`, and `control_group_resolution.json`. Every `crt_metadata.json` found is annotated with the resolved pool. |
| `sceptre_control_group.json` (optional) | Resolved `control_group`, `control_group_requested`, provenance, MOI and non-targeting cell count, written by `bin/inference_sceptre.R`. |
| `control_batch_composition.tsv`, `control_batch_composition_note.txt` | Written by `evaluate_controls.py` into the plots directory: cells and control-only cells per batch level, each level's control share and its share of all control cells, plus a note naming the level the controls concentrate in and a WARNING when one level holds >90% of them while holding <50% of the screen's cells. |
| `pipeline_dashboard.tar.gz`, `pipeline_qc_metrics.json` | Published at the top level of `--outdir` by `createDashboard`. The dashboard no longer publishes `inference_mudata.h5mu`, and it deletes `additional_qc/global_analysis/global_analysis_results.tsv` before packaging. |
| `DEMO_MODE_WARNING.txt` | Top level of `--outdir`, plus the filtered samplesheet in `pipeline_info/`, plus a copy inside `benchmark_output/` when `tf_benchmark` runs in demo mode. |
| `pipeline_info/` | New: `nextflow.config` (the last config *file* Nextflow loaded, copied verbatim, not a resolved or merged config), plus `nextflow.log`, `pipeline_manifest.config` (git commit/branch/remote, command line, Nextflow version, run name, success flag) and `original_samplesheet.csv\|tsv`. The timestamped `params_*.json` in that directory is unchanged from v1. `samplesheet.valid.csv` is no longer described in `docs/output.md`. |
| `seqspeccheck/` | `seqSpecCheck` now has a `publishDir` (`<outdir>/pipeline_outputs/seqspeccheck`); in v1 it had only `debug true`. |
| `cdna.txt`, `nascent_index.txt` | New unconditional declared outputs of `downloadReference`, extracted from the IGVF tarball and consumed by `mappingscRNA` for the `nac` workflow. |
| Repo docs (not run output) | Added: `docs/index.html`, `parameters.html`, `parameters_nfcore_style.html`, `mudata_schema.md`, `mudata_schema_catalog.tsv`, `qc_outputs.md`, `qc_outputs_flat.json`, `qc_outputs_flat.tsv`, `filtering_workflow.{md,html}`, `outputs/per_element_output_field_annotations.xlsx`, `AGENTS.md`. Modified, not added: `docs/README.md`, `docs/usage.md`, `docs/output.md`. |

`inference_mudata.h5mu` from `mergedResults` and `inference_sceptre` is now `optional: true`, and with the defaults (`INFERENCE_SCEPTRE_WRITE_MUDATA=false`, `INFERENCE_PERTURBO_WRITE_MUDATA=false`, and `mergedResults` called with `write_mudata=false` in the `default` workflow) those intermediates are not produced.

### Containers

| Entry | main | dev |
| --- | --- | --- |
| `containers.perturbo` | `ghcr.io/pinellolab/perturbo:sha-f3dc8ca` | `ghcr.io/pinellolab/perturbo@sha256:8c5d5a0005185e4744c88b338914b8eb7509ec6b90995eb75dc6266b993c417f` (v2.0 rc9) |
| `containers.base` | `…/conda-docker:latest` | `…/conda-docker:sha-91bc741` |
| `containers.sceptre` | `sjiang9/sceptre-igvf:0.1` | `sjiang9/sceptre-igvf:0.2` |
| `containers.aria2` | — | `biasofpriene/aria2c` (new; `downloadReference` moved here from `base`) |
| `containers.bediting` | — | `bioinfolucas/crispr_mapping:latest_0236` (new; both base-editing mappers) |
| `containers.gmmdemux` | `ghcr.io/lucassilvaferreira/gmm_demux_docker_crisprpipeline:latest` | unchanged |
| `containers.cleanser` | `…/cleanser:1.2.1` | unchanged |

Only the perturbo entry is a digest. `base` is the one image built from this repository (`docker-images/conda-docker/`, the tree's only recipe directory); `sceptre`, `aria2`, `bediting`, `cleanser` and `gmmdemux` are third-party or externally built tags whose contents nothing here characterises.

Other container-scope changes:

- All but one of `nextflow.config`'s seventeen process-scope container assignments, and all eight in `nextflow_cc.config`, are closures (`container = { params.containers.base }`) rather than eagerly-evaluated bare references, so a `-c` override of `params.containers.*` reaches the process directive. Seven of v1's eight assignments were converted; the nine added in this release were written as closures. The exception is `withName: 'demultiplex' { container = params.containers.gmmdemux }` (`nextflow.config:550`), still bare; overriding that one entry with `-c` does not take effect.
- The global `process.container` default (present in v1 as a bare reference) is now `{ params.containers.base }`, and `tests/test_nextflow_container_defaults.py::test_every_process_has_a_base_container_fallback` asserts that spelling stays in the global block.
- `sceptre_chunk_prepare|sceptre_chunk_merge` moved from `containers.base` to `containers.sceptre`.
- The selector `guide_assignment_sceptre|inference_sceptre` was split into two; `inference_perturbo|inference_perturbo_trans` collapsed to `inference_perturbo`.
- New selectors: `mappingGuideBaseEditing|mappingGuideBaseEditingFlash` -> `containers.bediting`, and `mergeMudata|buildCatalogElement|buildCatalogGuide` -> `containers.base`.
- `singularity.pullTimeout = '90 min'` added: the ~3.6 GB PerTurbo image can exceed Nextflow's 20-minute SIF-conversion default on a network filesystem.

### Base image recipe (`docker-images/conda-docker/nextflow.yaml`)

| Change | Consequence |
| --- | --- |
| `pyarrow` added | Parquet result tables and the Polars streaming paths. |
| `muon` **removed** | `bin/create_mdata.py` switched `from muon import MuData` to `from mudata import MuData`; no file in `bin/` or `modules/` imports muon at HEAD. |
| `anndata` pinned `>=0.11,<0.12` | Was unpinned. scanpy and numpy remain unpinned, so pinning the image freezes one solve of them. |

Together these are why `containers.base` had to move off `:latest`, which predates them and is rebuilt only by a manual `workflow_dispatch` of `conda_docker.yml`. `aria2` was already in this recipe in v1 and is unchanged.

### Schema-vs-config disagreements

`nextflow_schema.json` is not enforced, and it was itself substantially rewritten in this release (377 insertions / 57 deletions), so several of these are new schema values rather than stale ones. Eleven of its defaults disagree with `nextflow.config`, resolving the schema through its `$defs` and taking `nextflow.config`'s last-wins values. Restoring the run-specific defaults before release removed three further disagreements (`DEMO_MODE`, `reverse_complement_guides`, `QC_barcode_filter`), which now match the schema:

| Parameter | schema | config |
| --- | --- | --- |
| `GUIDE_ASSIGNMENT_capture_method` | `'crop-seq'` (v1 schema: `'CROP-seq'`) | `'CROP-seq'` (case only) |
| `QC_min_genes_per_cell` | `800` (changed here; v1 schema: `500`) | `500` |
| `ENABLE_BENCHMARK` | `true` | `false` |
| `INFERENCE_PERTURBO_GUIDE_GPU` | `'1'` | `'0'` |
| `GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold` | `0.8` | `'default'` |
| `GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep` | `5` | `'default'` |
| `spacer_tag` | `'GAGTACATGGGG'` | `''` |
| `is_10x3v3` | `true` | `false` |
| `INFERENCE_predefined_pairs_to_test` | `null` | `'path/to/file.csv'` |
| `REFERENCE_gtf_download_path` | GENCODE v43, moved **backward** in this release's schema rewrite; v1's schema and config both said v46 | GENCODE v46 |
| `INFERENCE_input_mudata` | `null` | `''` (a difference in spelling of "unset" rather than in value) |

Four of the eleven are one fossil, not four independent regressions. `c16ebd3` (12 Jun) rewrote `nextflow_schema.json` from `nextflow.config` as it stood that day, and `6a65d64` (20 Jul) then moved the config and left the schema alone. Every value the schema still advertises for `is_10x3v3` (`true`), `spacer_tag` (`'GAGTACATGGGG'`), `REFERENCE_gtf_download_path` (v43) and `QC_min_genes_per_cell` (`800`) is what `c16ebd3:nextflow.config` read at those keys; the first three were put there by `abc7d31` (14 May) and reverted by `6a65d64`. `QC_barcode_filter` was a fifth until the release restored the config to `'knee2'`, which is what the schema had said all along. So the schema's GENCODE v43 is a snapshot of a config state that no longer exists rather than a deliberate downgrade.

`REFERENCE_restrict_genes_to_gtf` is absent from the schema entirely, as are `max_cpus`, `max_memory` and `validate_params`.

### Open questions

These could not be settled from this checkout and are recorded rather than asserted.

- v1's `nextflow.config` is brace-unbalanced (the `process {` block opened at line 194 is never closed) and fails to parse under Nextflow 25.09.0-beta5957 with "Unexpected input: `{`". v1's parameters were resolved by appending one `}` to a throwaway copy. With that brace at EOF, v1's `singularity`, `tower`, `trace`, `report` and `timeline` scopes resolve as *process directives* (`process.singularity.enabled`, `process.trace.file`, …) rather than top-level scopes; dev has them at top level, so dev's `singularity.enabled`/`autoMounts`/`runOptions = '--nv'` genuinely take effect. Which Nextflow version, if any, ran v1's config as committed, and therefore whether v1.0.0's released behaviour was "singularity not enabled from the root config" or "config rejected at launch", could not be established.
- Whether the "Auto-injected by Pipeline Configurator" params blocks belong in the tracked `nextflow.config` at all. They still sit after the hand-written block and therefore win; the run-specific values they carried have been restored, but the mechanism remains, so the next configurator run can reintroduce the same problem.
- When `is_BaseEditing` went from `true` back to `false`. `baeea4d` (2025-11-14) added it as `true`; the transition to `false` appears as no diff hunk in `nextflow.config`'s history on any branch (it landed through a merge resolution), the parameter is absent from the file entirely on many commits in between, and the last commit where it demonstrably reads `true` is `c4b94a9` (2026-01-14). Whether `true` affected any real run in that window cannot be established from the repository either.
- The size of the result difference between the base-editing mappers and `mappingGuide` on the same input. Nothing in the repository compares them, and the conventions differ (CRISPR-Correct ambiguity spreading at Hamming 6 with `revcomp_protospacer=True`, or FLASH at `--tolerance`, versus kite pseudoalignment).
- Whether SCEPTRE really declines any factor with fifteen or more levels. The claim is a prose comment in `bin/inference_covariates.py:27-29` citing `MAX_N_LEVELS_ALLOWED` in the SCEPTRE package; nothing in this checkout verifies it against `sjiang9/sceptre-igvf:0.2`.
- The contents and provenance of `bioinfolucas/crispr_mapping:latest_0236` (including which `crispr_ambiguous_mapping` version it carries) and of `ghcr.io/lucassilvaferreira/gmm_demux_docker_crisprpipeline:latest`. Both are mutable-looking tags not built from anything here.
- `downloadReference` declares `path "cdna.txt"` and `path "nascent_index.txt"` unconditionally, but only the `use_igvf_reference = true` branch creates them; the else branch runs just `kb ref -d …`. A run with `use_igvf_reference = false` should therefore fail on missing outputs. Still the case at HEAD; not confirmed by running the engine.
- Whether `assets/barcode_replacement.txt` (73,728 rows) and `assets/barcode_3_2_1_bc1_replacement.tsv` (~442k rows) are the correct or complete maps for the CC-Perturb-seq chemistry. Only their shape is observable.
- kb-python's behaviour under the flags added for CC-Perturb-seq: whether `--sum total` and `--workflow standard` are accepted or ignored on the standard workflow, and whether `-r <file>` is what makes kb write the `counts_unfiltered_modified` directory every downstream reader now expects. Neither is checkable without the container.
- `mappingGuide` still passes `-w ${barcode_file}` twice on the `kb count` line, a leftover of the March rebase duplication. Whether kb tolerates the repeat or takes the last one was not verified.
- Whether `bin/perturbo_inference.py`, `bin/perturbo_inference_chunked.py` and `bin/chunk_mudata.py` are kept deliberately as reference or simply not yet deleted.
- Intent behind roughly 100 of the 211 commits is not verifiable: many subjects are opaque ("new_config", "small modifications", "last_changes_aprantly_run", "hto_new_fix", "BASE_EDITING") and carry no body, issue reference or test. Everything attributed to them here is observable behaviour read from the code. Author and commit dates also diverge badly across the March 2026 commits (author dates in W11-W13, commit dates on 2026-03-20/27 and 2026-04-16), so those ISO weeks do not always say when a change reached this branch.
- A few working artifacts are committed in the tree: `.DS_Store` (modified across three commits), a `.nextflow.pid` added by `1c26db8`, and Excel lock files (`~$per_element_output_field_annotations.xlsx`) added and removed several times. None affect behaviour.
- The 71 non-merge commits authored in Jan-Feb 2026 on `main` itself are outside this delta, so nothing here covers them. They are ancestors of dev, i.e. already released.

## v1.0.0 - 2025-AUG-19

- Initial release of IGVF Perturb-seq Pipeline
- Core pipeline functionality
- Documentation and examples

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
