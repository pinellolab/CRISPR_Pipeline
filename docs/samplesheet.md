# Samplesheet

The pipeline reads a CSV or TSV with one row per modality/read-pair (and, when relevant, lane) within a `measurement_sets` ID. RNA and guide rows belonging to the same biological measurement must share that ID. A hash row is required when `ENABLE_DATA_HASHING=true`. Do not substitute the old nf-core three-column samplesheet.

```csv
R1_path,R2_path,file_modality,measurement_sets,sequencing_run,lane,seqspec,barcode_onlist,guide_design,barcode_hashtag_map,cov_day
/data/SET1_RNA_R1.fastq.gz,/data/SET1_RNA_R2.fastq.gz,scRNA,SET1,RUN1,1,/refs/rna_seqspec.yml,/refs/barcodes.txt,/refs/guide_metadata.tsv,,0
/data/SET1_GUIDE_R1.fastq.gz,/data/SET1_GUIDE_R2.fastq.gz,gRNA,SET1,RUN1,1,/refs/guide_seqspec.yml,/refs/barcodes.txt,/refs/guide_metadata.tsv,,0
```

This is a **format example**, not real downloadable data. Paths must be accessible to every execution node. `R2_path` may be empty only if the assay and SeqSpec really use one read. RNA and guide files are often different libraries; never pair them solely by row order. The current validator requires nonempty `barcode_onlist` and `guide_design` values on **every** row, even though each is consumed by a specific modality; repeat the same appropriate path for the measurement set.

| Column | Meaning |
|---|---|
| `R1_path`, `R2_path` | FASTQ paths for that modality/read pair. Preserve correct R1/R2 orientation and lane pairing. |
| `file_modality` | `scRNA`, `gRNA`, or `hash`. |
| `measurement_sets` | Stable join key for RNA/guide/optional HTO from the same measurement set; used for per-set QC and batch labeling. |
| `sequencing_run`, `lane` | Provenance and lane labels; keep distinct lanes rather than silently combining unrelated files. |
| `seqspec` | SeqSpec YAML for this modality. |
| `barcode_onlist` | RNA cell-barcode whitelist/onlist, matching the chemistry. |
| `guide_design` | Guide metadata TSV/CSV for the guide modality. |
| `barcode_hashtag_map` | HTO metadata for `hash` rows when hashing is enabled. |
| `cov_*` | Optional measurement-set covariates, for example `cov_day`, `cov_replicate`, `cov_pool`. They must be consistent and nonempty on all rows of the set; see [covariates](samplesheet_covariates.md). |

The current parser requires the core columns shown in the example, except `R2_path` and `barcode_hashtag_map` are not required for non-hashing runs. Additional unrelated columns are accepted but are not automatically modeling covariates. `DEMO_MODE=true` deterministically chooses one complete eligible measurement set, retaining all of its lanes; its output is a diagnostic only. Download the [current format example](downloadable-examples/samplesheet.csv) and replace its illustrative paths before running.

Before launch, check that each set has RNA and guide rows; that the HTO row and mapping exist if requested; that `seqspec` and metadata files resolve; and that barcode onlist, barcode replacement (if used), and guide orientation agree across modalities. The dashboard's mapping overlap QC is the first run-time check of those assumptions.

For an HTO assay, compare the [three-modality example](downloadable-examples/hashing_samplesheet.csv) and set `ENABLE_DATA_HASHING=true`. The current validator still expects core fields on every HTO row; supply the same onlist and guide design path used by the set, as well as `barcode_hashtag_map`.
