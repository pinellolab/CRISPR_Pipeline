# Samplesheet reference

The pipeline is started with a CSV samplesheet passed via `--input`.
Each row represents **one set of FASTQ files for one modality of one sample**.
Multiple rows may share the same `measurement_sets` value — that is how
multi-modality experiments (RNA + gRNA ± hashing) are linked together.

## Example

```csv
R1_path,R1_md5sum,R2_path,R2_md5sum,file_modality,file_set,measurement_sets,sequencing_run,lane,flowcell_id,index,seqspec,barcode_onlist,onlist_method,strand_specificity,guide_design,barcode_hashtag_map
/path/to/R1.fastq.gz,md5checksum,/path/to/R2.fastq.gz,md5checksum,scRNA,sample_A,A,1,1,FLOWCELL01,,/path/to/rna_seqspec.yaml,/path/to/barcodes.tsv,no combination,3 prime to 5 prime,/path/to/guide_metadata.tsv,
/path/to/R1.fastq.gz,md5checksum,/path/to/R2.fastq.gz,md5checksum,gRNA,sample_A,A,1,1,FLOWCELL01,,/path/to/guide_seqspec.yaml,/path/to/barcodes.tsv,no combination,3 prime to 5 prime,/path/to/guide_metadata.tsv,
```

---

## Column descriptions

| Column | Required | Type | Allowed values | Description |
|--------|----------|------|----------------|-------------|
| `R1_path` | ✅ | string (file path or URL) | — | <!-- TODO: fill in description --> |
| `R1_md5sum` | ❌ | string | — | <!-- TODO: fill in description --> |
| `R2_path` | ❌ | string (file path or URL) | — | <!-- TODO: fill in description --> |
| `R2_md5sum` | ❌ | string | — | <!-- TODO: fill in description --> |
| `file_modality` | ✅ | string | `scRNA`, `gRNA`, `hash` | <!-- TODO: fill in description --> |
| `file_set` | ❌ | string | — | <!-- TODO: fill in description --> |
| `measurement_sets` | ✅ | string | — | <!-- TODO: fill in description — note: used as the batch/covariate grouping key; all rows belonging to the same experiment share this value --> |
| `sequencing_run` | ❌ | string | — | <!-- TODO: fill in description --> |
| `lane` | ❌ | string / integer | — | <!-- TODO: fill in description --> |
| `flowcell_id` | ❌ | string | — | <!-- TODO: fill in description --> |
| `index` | ❌ | string | — | <!-- TODO: fill in description --> |
| `seqspec` | ✅ | string (file path or URL) | — | <!-- TODO: fill in description — note: modality-specific seqspec YAML; one per modality per measurement set --> |
| `barcode_onlist` | ✅ | string (file path or URL) | — | <!-- TODO: fill in description --> |
| `onlist_method` | ❌ | string | `no combination`, … | <!-- TODO: fill in description --> |
| `strand_specificity` | ❌ | string | `3 prime to 5 prime`, … | <!-- TODO: fill in description --> |
| `guide_design` | ❌* | string (file path or URL) | — | <!-- TODO: fill in description — *required when `file_modality` is `gRNA` --> |
| `barcode_hashtag_map` | ❌* | string (file path or URL) | — | <!-- TODO: fill in description — *required when `file_modality` is `hash` --> |

---

## Channel fanout

The pipeline reads the samplesheet once with nf-schema and fans the validated
rows into per-modality channels immediately in
`subworkflows/local/utils_nfcore_perturbseq_pipeline/main.nf`:

| Channel | Filter |
|---------|--------|
| `ch_rna` | `file_modality == 'scrna'` |
| `ch_guide` | `file_modality == 'grna'` |
| `ch_hash` | `file_modality == 'hash'` |

Per-modality auxiliary files (`seqspec`, `barcode_onlist`, `guide_design`,
`barcode_hashtag_map`) are emitted as dedicated single-value channels so
downstream subworkflows receive `path` inputs rather than having to unpack
them from the meta map themselves.
