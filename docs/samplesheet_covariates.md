# Samplesheet and checkpoint covariates

These options are included in the current `dev` pipeline. The two paths serve different scopes: `cov_*` columns describe an entire measurement set; `external_*` columns introduced at a checkpoint describe individual retained cells.

Add columns whose names begin with `cov_` to the usual CSV/TSV samplesheet. For example:

```csv
R1_path,R2_path,file_modality,measurement_sets,sequencing_run,lane,seqspec,barcode_onlist,guide_design,barcode_hashtag_map,cov_day,cov_replicate,cov_pool
/reads/A_R1.fastq.gz,/reads/A_R2.fastq.gz,scRNA,SET_A,RUN_A,1,/ref/rna.yaml,/ref/barcodes.txt,/ref/guides.tsv,,0,rep1,poolA
/reads/A_R1.fastq.gz,/reads/A_R2.fastq.gz,gRNA,SET_A,RUN_A,1,/ref/guide.yaml,/ref/barcodes.txt,/ref/guides.tsv,,0,rep1,poolA
/reads/B_R1.fastq.gz,/reads/B_R2.fastq.gz,scRNA,SET_B,RUN_B,1,/ref/rna.yaml,/ref/barcodes.txt,/ref/guides.tsv,,7,rep2,poolB
/reads/B_R1.fastq.gz,/reads/B_R2.fastq.gz,gRNA,SET_B,RUN_B,1,/ref/guide.yaml,/ref/barcodes.txt,/ref/guides.tsv,,7,rep2,poolB
```

Every row of one `measurement_sets` ID must give the **same, nonempty** value for each declared `cov_*` column, including across RNA and guide modalities. Other extra columns are allowed but not used as covariates. Numeric columns such as `cov_day=0,7` are treated as continuous; labels such as `day0,day7`, `rep1,rep2`, and `poolA,poolB` are categorical. Labels with leading zeros remain categorical. This path represents **measurement-set-level** metadata; for per-cell covariates or external cell filtering, use the [pre-inference checkpoint continuation](preinference_checkpoint_user_guide.md).

`cov_*` values are copied to RNA cell metadata and appear as PCA/UMAP color overlays. They do not change PCA coordinates or normalized expression. With `INFERENCE_COV_CUSTOM = 'auto'` (default), the pipeline tries to include every `cov_*` column in both SCEPTRE and PerTurbo. Use `'none'` to disable all in inference while retaining PCA overlays, or a comma-separated list such as `'cov_day,cov_pool'`. Categorical columns are reference-coded into numeric indicators to support PerTurbo's one-batch-factor CLI; the same indicators enter SCEPTRE. Missing, constant, high-cardinality, or collinear requested covariates cause a clear pre-inference failure.

Important: `batch` is the measurement set. A `cov_day`, `cov_replicate`, or `cov_pool` value assigned once per measurement set is therefore **collinear with full batch adjustment**. To estimate a set-level covariate, use a scientifically justified configuration with `INFERENCE_COV_BATCH = false`, or redesign the batch variable. The pipeline must not silently drop the custom variable from only one method. Turning batch off changes the statistical question; compare such runs as sensitivity analyses, not as a direct method-only comparison.

Available shared inference switches (all named `INFERENCE_COV_*`):

| Switch | Default | Quantity |
| --- | --- | --- |
| `MITOCHONDRIAL` | `false` | Optional pre-feature-filter mitochondrial percentage |
| `BATCH` | `true` | Measurement-set factor |
| `RIBOSOMAL` | `false` | Pre-feature-filter ribosomal percentage; must exist in cell metadata |
| `DETECTED_GENES` | `true` | Log detected RNA genes from the full matrix; PR #101 default |
| `GUIDE_COUNT` | `false` | Number of assigned guides from `guide_assignment` |
| `GUIDE_UMIS` | `true` | Log total raw guide UMIs; PR #101 default |
| `CUSTOM` | `'auto'` | All, none, or selected `cov_*` columns |

The PR #101 default design is log guide UMI depth, log RNA UMI depth, log detected genes, and batch. Disabling `DETECTED_GENES` or `GUIDE_UMIS` is a sensitivity analysis. Guide count is opt-in. RNA library size remains a model depth/offset term and is not disabled by these switches. The original [SCEPTRE paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02545-2) considered RNA depth, detected genes, detected gRNAs, mitochondrial percentage, and batch; its [formula guidance](https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#formula) notes that default covariates differ by MOI. Ribosomal fraction is an additional pipeline option, not a SCEPTRE-paper default.

At the pre-inference checkpoint, per-cell numeric columns named `external_*` may be added to `gene.obs` and listed in `mdata.uns['external_covariates']`. On continuation they are validated, copied to top-level cell metadata, and included in both inference methods. The cell barcode index must match the checkpoint exactly; this does not change guide assignment or the count matrix. An `external_*` term must be finite and nonconstant. The pipeline checks the combined design for collinearity before fitting.

For a runnable checkpoint example, including the TSV keyed by full cell IDs and the validation rules, see [pause, edit, continue](preinference_checkpoint_user_guide.md). The same covariate design is passed to both methods; always review the resolved inference config and report the model design with the result.
