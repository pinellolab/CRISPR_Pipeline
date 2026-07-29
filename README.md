![Alt text](https://github.com/pinellolab/CRISPR_Pipeline/blob/main/images/crispr_pipeline.png)


[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/crispr)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23crispr-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/crispr)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

# CRISPR Pipeline

A comprehensive pipeline for single-cell Perturb-Seq analysis that enables robust processing and analysis of CRISPR screening data at single-cell resolution.


## Documentation Links

- [Input parameters and explanations](https://pinellolab.github.io/CRISPR_Pipeline/)
- [Colab seqspec checker](https://colab.research.google.com/drive/1IfSU9Oxf3-JIuTt8oGfWjCqsiv6ff4LB#scrollTo=XMxGDEuBo1w9): Check whether your seqspec is extracting cell barcodes, UMIs, guides, hashes, and transcripts correctly.

## Prerequisites

Nextflow and Singularity must be installed before running the pipeline:

### Nextflow (version > 24)
Workflow manager for executing the pipeline:

```bash
conda install bioconda::nextflow
```
### Singularity
Container platform that must be available in your execution environment.

### Nextflow Tower Integration
This is a seamless pipeline execution monitoring system that offers a web-based interface for workflow management.

To enable Nextflow Tower, we require a TOWER_ACCESS_TOKEN.

To obtain your token:
1. Create/login to your account at cloud.tower.nf
2. Navigate to Settings > Your tokens
3. Click "Add token" and generate a new token
4. Set as environment variable: `export TOWER_ACCESS_TOKEN=your_token_here`

## Pipeline Installation

To install the pipeline:

```bash
git clone https://github.com/pinellolab/CRISPR_Pipeline.git
```

## Input Requirements

### File Descriptions

#### FASTQ Files
- `{sample}_R1.fastq.gz`: Contains cell barcode and UMI sequences
- `{sample}_R2.fastq.gz`: Contains transcript sequences

#### YAML Configuration Files (see example_data/)
- `rna_seqspec.yml`: Defines RNA sequencing structure and parameters
- `guide_seqspec.yml`: Specifies guide RNA detection parameters
- `hash_seqspec.yml`: Defines cell hashing structure (required if using cell hashing)
- `barcode_onlist.txt`: List of valid cell barcodes

#### Metadata Files (see example_data/)
- `guide_metadata.tsv`: Contains guide RNA information and annotations
- `hash_metadata.tsv`: Cell hashing sample information (required if using cell hashing)

For detailed specifications, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

## Running the Pipeline 

### Pipeline Configuration

Before running the pipeline, customize the configuration files for your environment:

#### 1. Data and Analysis Parameters (`nextflow.config`)

Update the pipeline-specific parameters in the `params` section. The tables below document the options that change biological processing, guide assignment, inference, filtering, or optional analysis outputs. Compute resources, containers, cloud settings, email hooks, and other execution-only settings are documented separately in the compute configuration section.

The same fields are also represented in `nextflow_schema.json`, using nf-core-style parameter groups so GUI/launch tools can show the options with descriptions.

Runtime/debug/internal keys such as `DEBUG_VAR`, dashboard asset paths (`css`, `js`, `svg`), and currently unused placeholders such as `INFERENCE_SCEPTRE_formula_object` are intentionally not listed here.

##### Input and output options

| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `input` | `null` | CSV or TSV samplesheet path | Samplesheet containing modality-specific FASTQ paths and metadata used to create the RNA, guide, and optional hashing channels. |
| `outdir` | `./pipeline_outputs` | Directory path | Output directory for published final files, including MuData, per-guide/per-element result tables, dashboard archive, and QC metrics. |
| `DEMO_MODE` | `false` | `true`, `false` | Pre-run only. Selects the complete `measurement_sets` group with the fewest scRNA FASTQ files. If hashing is enabled, hash rows are required and retained too. Demo outputs are not final results. |

##### Demo pre-runs

Use demo mode to validate configuration, references, SeqSpecs, and pipeline execution on one matched measurement set before starting the full run:

```bash
nextflow run . \
  --input samplesheet.csv \
  --outdir demo_outputs \
  --DEMO_MODE true \
  -profile local
```

Selection is deterministic across repeated runs and independent of samplesheet row order. The pipeline finds every `measurement_sets` ID containing both `scRNA` and `gRNA`, counts its distinct non-empty scRNA `R1_path` and `R2_path` values, and selects the set with the fewest scRNA FASTQ files. Ties are resolved using the `measurement_sets` ID. It keeps every supported RNA/guide/hash row for that ID, including multiple lanes/read pairs. When `ENABLE_DATA_HASHING = true`, only sets that also contain `hash` are eligible. The filtered file is saved as `pipeline_info/demo_samplesheet.csv` or `.tsv`; the selected ID and scRNA FASTQ count are recorded in `DEMO_MODE_WARNING.txt`.

Demo mode is strictly for pre-runs. The pipeline prints warnings at startup and completion, writes `DEMO_MODE_WARNING.txt` at the output root and in `pipeline_info/`, marks the dashboard benchmark card, and adds a red pre-run warning to the TF benchmark plot plus its benchmark output directory. Rerun with `DEMO_MODE = false` for final results.

##### Assay and library options

[input parameters and explanations] (https://pinellolab.github.io/CRISPR_Pipeline/) 



| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `ENABLE_DATA_HASHING` | `false` | `true`, `false` | Enables the hashing workflow: hash seqspec checks, hash mapping, hashtag filtering, demultiplexing, hash-aware MuData creation, and hash dashboard sections. |
| `ENABLE_SCRUBLET` | `false` | `true`, `false` | Runs Scrublet doublet detection before guide assignment in the non-hashing workflow. |
| `is_10x3v3` | `true` | `true`, `false` | Controls 10x Genomics 3' v3 feature-barcode chemistry (`10XV3`, `kite:10xFB`) for guide or hashing mapping depending on `ENABLE_DATA_HASHING`. RNA mapping always uses the RNA seqspec. Case 1: when `ENABLE_DATA_HASHING = false` and `is_10x3v3 = true`, guide mapping uses the 10x v3 feature-barcode kb settings instead of deriving guide chemistry from the guide seqspec. Case 2: when `ENABLE_DATA_HASHING = true` and `is_10x3v3 = true`, guide and RNA mapping use their seqspecs, while hash/HTO mapping uses the 10x v3 feature-barcode kb settings. This second case supports 10x v3 HTO data where barcode replacement/translation may be needed so hash, RNA, and guide barcodes match downstream. |
| `reverse_complement_guides` | `false` | `true`, `false` | Reverse-complements guide spacer sequences while building the guide reference, preserving the metadata fields. |
| `spacer_tag` | `GAGTACATGGGG` | DNA sequence, empty string, or `null` | Recommended 12 bp sequence immediately upstream of the guide spacer. When provided, guide mapping searches the whole guide read around this tag instead of relying only on fixed seqspec feature coordinates. |
| `scrna_workflow` | `standard` | `standard`, `nac` | Selects the kb count RNA workflow. `standard` performs mature transcript counting; `nac` performs nascent-aware counting with cDNA and nascent references for unspliced/nascent signal. |
| `use_multimapping` | `false` | `true`, `false` | Passes kb count multimapping mode for scRNA mapping and keeps that setting during AnnData concatenation. |
| `replace_barcodes` | `false` | `true`, `false` | Enables the CC-Perturb-seq barcode replacement strategy during RNA and guide mapping. When enabled, kb receives the replacement table and downstream concatenation reads `counts_unfiltered_modified`. |
| `bc_replacement_file` | `''` | File path | Replacement table used by kb count when `replace_barcodes = true`. |
| `DUAL_GUIDE` | `false` | `true`, `false` | Enables dual-guide-aware aggregation when concatenating guide-assigned MuData, for dual-guide perturbation designs. |
| `Multiplicity_of_infection` | `high` | `high`, `low` | Records screen MOI in MuData. `high` allows multi-guide cell contexts; `low` records a mostly zero-or-one-guide design and changes how guide assignments and perturbation tests should be interpreted. |

##### Reference options

| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `use_igvf_reference` | `true` | `true`, `false` | Uses the prebuilt IGVF transcriptome reference for RNA mapping. When enabled, `REFERENCE_transcriptome` is ignored for RNA reference selection. |
| `REFERENCE_transcriptome` | `human` | kb-python transcriptome name | Transcriptome reference name used by kb-python only when `use_igvf_reference = false`; ignored for RNA mapping when the IGVF reference is enabled. |
| `REFERENCE_gtf_download_path` | GENCODE v43 URL | URL | GTF annotation URL used only when no local GTF exists. Reference-selection settings are ignored when `use_igvf_reference = true`, but the GTF is still needed for gene annotation and cis pair construction. |
| `REFERENCE_gtf_local_path` | `/path/to/gencode_gtf.gtf.gz` | File path | Local GTF annotation. If present, preprocessing/inference use it instead of downloading `REFERENCE_gtf_download_path`; RNA reference selection still follows `use_igvf_reference`. |

##### Quality control options

The complete machine-readable QC output catalog is available as
[flat JSON](docs/qc_outputs_flat.json), [TSV](docs/qc_outputs_flat.tsv), and a
[rendered table](docs/qc_outputs.md). Regenerate all three with
`bin/export_qc_metric_catalog.py`.

| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `QC_min_genes_per_cell` | `800` | Integer | Minimum detected genes required to keep a cell when `QC_barcode_filter = 'none'`. A gene is counted as present in a cell when its RNA count is greater than zero. |
| `QC_min_cells_per_gene` | `0.05` | Fraction in `[0, 1)` | Minimum retained-cell fraction required to keep a gene during guide-assignment aggregation. `0` retains every gene detected in at least one cell. |
| `TAPSEQ_QC_MODE` | `false` | `true`, `false` | TAP-seq gene-retention mode. It removes the standard 10-cell preprocessing floor, retaining every observed gene before the final fractional support filter. Use a small fraction such as `0.000001` when all observed TAP-seq genes should be retained. |
| `QC_pct_mito` | `15` | `0` to `100` | Maximum mitochondrial read percentage allowed per cell during preprocessing. |
| `QC_batch_col` | `batch` | Observation column name | Batch column used in additional QC plots. |
| `QC_barcode_filter` | `knee2` | `none`, `knee`, `knee2` | RNA barcode filtering strategy based on total cell RNA UMIs. `knee` uses the first barcode-rank knee and is more permissive; `knee2` searches the high-UMI segment before knee1 for a second, stricter knee; `none` skips UMI-knee filtering and applies `QC_min_genes_per_cell`. If the requested knee cannot be found, barcode filtering is skipped and the min-gene filter is not applied. |

##### Guide assignment options

| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `GUIDE_ASSIGNMENT_method` | `sceptre` | `sceptre`, `cleanser` | Selects the guide-to-cell assignment method before inference. |
| `GUIDE_ASSIGNMENT_capture_method` | `crop-seq` | `crop-seq`, `direct-capture` | Recorded in MuData and passed directly to CLEANSER as `--crop-seq` or `--direct-capture` when `GUIDE_ASSIGNMENT_method = 'cleanser'`. |
| `GUIDE_ASSIGNMENT_cleanser_probability_threshold` | `1` | `0` to `1` | Probability threshold used by Cleanser guide assignment. |
| `GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold` | `0.8` | `0` to `1` | Posterior probability threshold for SCEPTRE mixture-based guide assignment. |
| `GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep` | `5` | Integer `>= 1` | Number of EM initializations used by SCEPTRE guide assignment. |

##### Inference options

| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `INFERENCE_method` | `default` | `default`, `sceptre`, `perturbo`, `sceptre,perturbo` | Selects inference workflow. `default` runs cis SCEPTRE, cis PerTurbo, trans PerTurbo, and writes merged cis/trans outputs. |
| `INFERENCE_input_mudata` | `null` | MuData `.h5mu` path | Required only for `-entry INFERENCE_FROM_MUDATA`, where inference is rerun from an existing post-guide-assignment MuData file. |
| `INFERENCE_target_guide_pairing_strategy` | `default` | `default`, `by_distance`, `predefined_pairs` | Controls how guide-target pairs are built before inference. `default` builds the standard cis pairs; `by_distance` uses genomic distance; `predefined_pairs` uses a user-provided table. |
| `INFERENCE_predefined_pairs_to_test` | `null` | CSV path | Pair table required when `INFERENCE_target_guide_pairing_strategy = 'predefined_pairs'`. |
| `INFERENCE_max_target_distance_bp` | `1000000` | Integer bp distance | Maximum guide-target genomic distance used for cis pair construction. |
| `INFERENCE_PERTURBO_BATCH_SIZE` | `4096` | Integer `>= 1` | Batch size passed to PerTurbo inference in cis and trans runs. |
| `INFERENCE_PERTURBO_TRANS_MAX_GENES_PER_CHUNK` | `8000` | Integer; `<= 0` disables chunking | Maximum genes per chunk for trans PerTurbo all-by-all inference. |
| `INFERENCE_SCEPTRE_side` | `both` | `both`, `left`, `right` | Alternative-hypothesis side passed to SCEPTRE inference. |
| `INFERENCE_SCEPTRE_grna_integration_strategy` | `union` | SCEPTRE strategy string | Guide RNA integration strategy passed to SCEPTRE inference. |
| `INFERENCE_SCEPTRE_resampling_approximation` | `skew_normal` | SCEPTRE approximation string | Resampling approximation passed to SCEPTRE inference. |
| `INFERENCE_SCEPTRE_control_group` | `complement` | SCEPTRE control group string | Control group strategy passed to SCEPTRE inference. |
| `INFERENCE_SCEPTRE_resampling_mechanism` | `default` | SCEPTRE mechanism string | Resampling mechanism passed to SCEPTRE inference. |
| `INFERENCE_SCEPTRE_CHUNK_MODE` | `auto` | `auto`, `off`, `force` | Controls gene chunking before SCEPTRE inference. `auto` chunks large matrices, `off` keeps a single input, and `force` chunks regardless of matrix size. |
| `INFERENCE_SCEPTRE_MAX_MATRIX_ENTRIES` | `2147483647` | Integer `>= 1` | Cell-by-gene matrix size threshold used by SCEPTRE auto chunking. |
| `INFERENCE_SCEPTRE_GENE_CHUNK_SIZE` | `1000` | Integer `>= 1` | Number of genes per SCEPTRE chunk when chunking is enabled. |
| `INFERENCE_SCEPTRE_FORCE_CHUNK` | `false` | `true`, `false` | Compatibility flag that passes `--force-chunk` to the SCEPTRE chunking step. |

##### Optional analysis and display options

| Parameter | Default | Options | Pipeline context |
|---|---:|---|---|
| `ENABLE_BENCHMARK` | `true` | `true`, `false` | Runs the optional transcription-factor benchmark workflow after inference. In demo mode, benchmark plots, files, and the dashboard block are marked as pre-run-only and not final results. |
| `ENCODE_BED_DIR` | `${projectDir}/encode_bed_files` | Directory path | Directory containing ENCODE BED files used by the optional benchmark workflow. |
| `NETWORK_custom_central_nodes` | `undefined` | Comma-separated node names or `undefined` | Custom central nodes to highlight in network plots. |
| `NETWORK_central_nodes_num` | `1` | Integer `>= 0` | Number of central nodes to highlight when custom central nodes are not provided. |

### Run default inference from an existing MuData

Use this mode when you already have a post-guide-assignment MuData file (for example, after custom cell/guide filtering) and only need to rerun the default inference workflow.

```bash
nextflow run main.nf \
  -entry INFERENCE_FROM_MUDATA \
  -profile local \
  --INFERENCE_input_mudata /path/to/filtered_input.h5mu \
  --INFERENCE_method default \
  --INFERENCE_target_guide_pairing_strategy default \
  --REFERENCE_gtf_local_path /path/to/gencode_gtf.gtf.gz \
  --outdir ./outputs_mudata_inference
```

Notes:
- This entrypoint runs inference only (no mapping, guide assignment, dashboard, or additional QC workflows).
- GTF is still required in default mode because the pipeline constructs cis pairs before running inference.
- If `REFERENCE_gtf_local_path` does not exist, the pipeline uses `REFERENCE_gtf_download_path`.

#### 2. Compute Environment Configuration

Choose and configure your compute profile by updating the relevant sections:

##### 🖥️ **Local**
```groovy
// Resource limits (adjust based on your machine)
max_cpus = 8           // Number of CPU cores available
max_memory = '32.GB'   // RAM available for the pipeline

// Run with: nextflow run main.nf -profile local
```

##### 🏢 **SLURM Cluster**
```groovy
// Resource limits (adjust based on cluster specs)
max_cpus = 128
max_memory = '512.GB'

// Update SLURM partitions in profiles section:
slurm {
    process {
        queue = 'short,normal,long'  // Replace with your partition names
    }
}

// Run with: nextflow run main.nf -profile slurm
```

##### ☁️ **Google Cloud Platform**
```groovy
// Update GCP settings
google_bucket = 'gs://your-bucket-name'
google_project = 'your-gcp-project-id'
google_region = 'us-central1'  // Choose your preferred region

// Resource limits
max_cpus = 128
max_memory = '512.GB'

// Run with (see more in GCP_user_notebook.ipynb): 
// export GOOGLE_APPLICATION_CREDENTIALS="/path/to/your/pipeline-service-key.json"
// nextflow run main.nf -profile google 
```

#### 3. Container Configuration

The pipeline uses pre-built containers. Update if you have custom versions:

```groovy
containers {
   base     = 'ghcr.io/pinellolab/crispr_pipeline/conda-docker'
   cleanser = 'ghcr.io/gersbachlab-bioinformatics/cleanser:1.2.1'
   sceptre  = 'sjiang9/sceptre-igvf:0.1'
   perturbo = 'ghcr.io/pinellolab/perturbo'
}
```

## 🎯 Resource Sizing Guidelines

### Recommended Starting Values:

| Environment | max_cpus | max_memory | Notes |
|------------|----------|------------|--------|
| **Local (development)** | 4-8 | 16-32GB | For testing small datasets |
| **Local (full analysis)** | 8-16 | 64-128GB | For complete runs |
| **SLURM cluster** | 64-128 | 256-512GB | Adjust based on node specs |
| **Google Cloud** | 128+ | 512GB+ | Can scale dynamically |

## 🔧 Testing Your Configuration

1. **Validate syntax:**
   ```bash
   nextflow config -profile local  # Test local profile
   nextflow config -profile slurm  # Test SLURM profile
   ```

2. **Test with small dataset:**
   ```bash
   # Start with a subset of your data
   # Make all scripts executable (required for pipeline execution)
   chmod +x bin/*
   # RUN THE PIPELINE
   nextflow run main.nf -profile local --input small_test.tsv -outdir ./Outputs
   ```


## 💡 Pro Tips

- **Start conservative:** Begin with lower resource limits and increase as needed
- **Profile-specific limits:** The pipeline automatically scales resources based on retry attempts
- **Development workflow:** Use local profile for code testing, cluster/cloud for production runs

## 🚨 Common Issues

- **Memory errors:** Increase `max_memory` if you see out-of-memory failures
- **Queue timeouts:** Adjust SLURM partition names to match your cluster
- **Permission errors:** Ensure your Google Cloud service account has proper permissions
- **Container issues:** Verify Singularity is available on your system
- **Missing files**: Double-check paths in `nextflow.config` and actual files in `example_data`

## Output Description

All paths below are relative to the directory supplied with `--outdir`.

Final inference artifacts are written once, under `pipeline_outputs/`. The dashboard directory is visualization-only and does not contain duplicate copies of `inference_mudata.h5mu` or the cis/trans TSV outputs.

### Final inference outputs

Within `pipeline_outputs/`, you will find:

| File | Description |
|---|---|
| `inference_mudata.h5mu` | Final MuData object containing processed modalities and inference results. |
| `cis_per_element_output.tsv.gz` | Cis element-level inference results. |
| `cis_per_guide_output.tsv.gz` | Cis guide-level inference results. |
| `trans_per_element_output.tsv.gz` | Trans element-level inference results. |
| `trans_per_guide_output.tsv.gz` | Trans guide-level inference results. |
| `catalog_per_element_output.tsv.gz` | Catalog-formatted per-element table merging cis SCEPTRE and trans PerTurbo results. |

All result tables are tab-separated and gzip-compressed.

### Cis analysis

The cis outputs report guide-gene or target-element-gene tests restricted to the configured cis pairing strategy.

| File | Description |
|---|---|
| `cis_per_guide_output.tsv.gz` | Inference results for cis guide-gene pairs with guides tested independently. |
| `cis_per_element_output.tsv.gz` | Inference results for cis element-gene pairs with guides grouped by intended target fields. |

#### `cis_per_guide_output.tsv.gz`

| Column | Description |
|---|---|
| `gene_id` | ENSEMBL gene ID |
| `guide_id` | Guide identifier (guide name) |
| `sceptre_log2_fc` | SCEPTRE effect size estimate (log2 fold-change) |
| `sceptre_p_value` | SCEPTRE (uncorrected) p_value  of differential expression |
| `sceptre_q_value` | BH-adjusted SCEPTRE p-value. |
| `sceptre_fc_se` | SCEPTRE fold-change standard error. |
| `sceptre_negLog10p` | SCEPTRE significance score: `-log10(max(sceptre_p_value, 1e-300))`. |
| `sceptre_log10_p_value` | Backward-compatible alias of `sceptre_negLog10p`. |
| `perturbo_log2_fc` | PerTurbo effect size estimate (log2 fold-change) |
| `perturbo_p_value` | PerTurbo (uncorrected) posterior probability  of differential expression |
| `perturbo_q_value` | BH-adjusted PerTurbo p-value, computed within this cis output table. |
| `perturbo_fc_se` | PerTurbo posterior standard error for the log2 fold-change estimate. |
| `perturbo_negLog10p` | PerTurbo significance score: `-log10(max(perturbo_p_value, 1e-300))`. |
| `perturbo_log10_p_value` | Backward-compatible alias of `perturbo_negLog10p`. |

#### `cis_per_element_output.tsv.gz`

| Column | Description |
|---|---|
| `gene_id` | ENSEMBL gene ID |
| `intended_target_name` | Intended target element name |
| `intended_target_chr` | Intended target chromosome |
| `intended_target_start` | Intended target start coordinate |
| `intended_target_end` | Intended target end coordinate |
| `sceptre_log2_fc` | SCEPTRE effect size estimate (log2 fold-change) |
| `sceptre_p_value` | SCEPTRE (uncorrected) p_value |
| `sceptre_q_value` | BH-adjusted SCEPTRE p-value. |
| `sceptre_fc_se` | SCEPTRE fold-change standard error. |
| `sceptre_negLog10p` | SCEPTRE significance score: `-log10(max(sceptre_p_value, 1e-300))`. |
| `sceptre_log10_p_value` | Backward-compatible alias of `sceptre_negLog10p`. |
| `perturbo_log2_fc` | PerTurbo effect size estimate (log2 fold-change) |
| `perturbo_p_value` | PerTurbo (uncorrected) posterior probability  of differential expression |
| `perturbo_q_value` | BH-adjusted PerTurbo p-value, computed within this cis output table. |
| `perturbo_fc_se` | PerTurbo posterior standard error for the log2 fold-change estimate. |
| `perturbo_negLog10p` | PerTurbo significance score: `-log10(max(perturbo_p_value, 1e-300))`. |
| `perturbo_log10_p_value` | Backward-compatible alias of `perturbo_negLog10p`. |
| `element_id` | Element identifier (equal to `element_name` in this pipeline). |
| `element_type` | Element type derived from guide metadata (`guide.var['type']`). |
| `element_chr` | Element chromosome. |
| `element_start` | Element start coordinate. |
| `element_end` | Element end coordinate. |
| `element_name` | Element name mapped from `intended_target_name`. |
| `guide_ids` | Sorted unique guide IDs for the element, separated by `;`. |
| `gene_name` | Gene symbol from local gene metadata when available. |
| `nPerturbedCells` | Number of unique cells assigned at least one guide for the element. |

`intended_target_name` for non-targeting controls is bucketed as `non-targeting|N` (for example, `non-targeting|1`).
SCEPTRE outputs contain discovery-analysis results only (calibration-check rows are not exported).

### Trans analysis

The trans outputs report PerTurbo all-by-all trans tests.

| File | Description |
|---|---|
| `trans_per_guide_output.tsv.gz` | PerTurbo inference results for all guide-gene pairs, with guides tested independently. |
| `trans_per_element_output.tsv.gz` | PerTurbo inference results for all element-gene pairs, grouped by intended target fields. |

#### `trans_per_guide_output.tsv.gz`

| Column | Description |
|---|---|
| `gene_id` | ENSEMBL gene ID |
| `guide_id` | Guide identifier (guide name). |
| `perturbo_log2_fc` | PerTurbo effect size (log2 fold-change) |
| `perturbo_p_value` | PerTurbo (uncorrected) posterior probability of differential expression |
| `perturbo_q_value` | BH-adjusted PerTurbo p-value, computed within this trans output table. |
| `perturbo_fc_se` | PerTurbo posterior standard error for the log2 fold-change estimate. |
| `perturbo_negLog10p` | PerTurbo significance score: `-log10(max(perturbo_p_value, 1e-300))`. |
| `perturbo_log10_p_value` | Backward-compatible alias of `perturbo_negLog10p`. |

#### `trans_per_element_output.tsv.gz`

| Column | Description |
|---|---|
| `gene_id` | ENSEMBL gene ID |
| `intended_target_name` | Intended target element name. |
| `intended_target_chr` | Intended target chromosome. |
| `intended_target_start` | Intended target start coordinate. |
| `intended_target_end` | Intended target end coordinate. |
| `perturbo_log2_fc` | PerTurbo effect size (log2 fold-change) |
| `perturbo_p_value` | PerTurbo (uncorrected) posterior probability of differential expression |
| `perturbo_q_value` | BH-adjusted PerTurbo p-value, computed within this trans output table. |
| `perturbo_fc_se` | PerTurbo posterior standard error for the log2 fold-change estimate. |
| `perturbo_negLog10p` | PerTurbo significance score: `-log10(max(perturbo_p_value, 1e-300))`. |
| `perturbo_log10_p_value` | Backward-compatible alias of `perturbo_negLog10p`. |
| `element_id` | Element identifier (equal to `element_name` in this pipeline). |
| `element_type` | Element type derived from guide metadata (`guide.var['type']`). |
| `element_chr` | Element chromosome. |
| `element_start` | Element start coordinate. |
| `element_end` | Element end coordinate. |
| `element_name` | Element name mapped from `intended_target_name`. |
| `guide_ids` | Sorted unique guide IDs for the element, separated by `;`. |
| `gene_name` | Gene symbol from local gene metadata when available. |
| `nPerturbedCells` | Number of unique cells assigned at least one guide for the element. |

#### `catalog_per_element_output.tsv.gz`

| Column | Description |
|---|---|
| `sceptre_log2_fc` | SCEPTRE effect size estimate from cis per-element results. |
| `sceptre_p_value` | SCEPTRE p-value from cis per-element results. |
| `sceptre_q_value` | BH-adjusted SCEPTRE p-value from cis per-element results. |
| `sceptre_fc_se` | SCEPTRE fold-change standard error from cis per-element results. |
| `sceptre_negLog10p` | Catalog-facing SCEPTRE significance score: `-log10(max(sceptre_p_value, 1e-300))`. Prefer this for catalog/export consumers that should avoid extremely small raw p-values. |
| `sceptre_log10_p_value` | Backward-compatible alias of `sceptre_negLog10p`. |
| `perturbo_log2_fc` | PerTurbo effect size estimate from trans per-element results. |
| `perturbo_p_value` | PerTurbo p-value from trans per-element results. |
| `perturbo_q_value` | BH-adjusted PerTurbo p-value from trans per-element results. |
| `perturbo_fc_se` | PerTurbo posterior standard error from trans per-element results. |
| `perturbo_negLog10p` | Catalog-facing PerTurbo significance score: `-log10(max(perturbo_p_value, 1e-300))`. Prefer this for catalog/export consumers that should avoid extremely small raw p-values. |
| `perturbo_log10_p_value` | Backward-compatible alias of `perturbo_negLog10p`. |
| `element_id` | Element identifier (equal to `element_name` in this pipeline). |
| `element_type` | Element type derived from guide metadata (`guide.var['type']`). |
| `element_chr` | Element chromosome. |
| `element_start` | Element start coordinate. |
| `element_end` | Element end coordinate. |
| `element_name` | Element name mapped from `intended_target_name`. |
| `guide_ids` | Sorted unique guide IDs for the element, separated by `;`. |
| `gene_name` | Gene symbol from local gene metadata when available. |
| `gene_id` | ENSEMBL gene ID. |
| `nPerturbedCells` | Number of unique cells assigned at least one guide for the element. |

For details, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

### Pipeline dashboard

Within `pipeline_dashboard/`, you will find the interactive dashboard and supporting visualization files. A compressed copy of this directory is also written to the top level of `--outdir` as `pipeline_dashboard.tar.gz`.

The dashboard directory and archive intentionally do not include `inference_mudata.h5mu`, `cis_per_element_output.tsv.gz`, `cis_per_guide_output.tsv.gz`, `trans_per_element_output.tsv.gz`, `trans_per_guide_output.tsv.gz`, or `catalog_per_element_output.tsv.gz`; use the copies in `pipeline_outputs/` as the single source of final analysis outputs.

The pipeline produces several figures:

1. **Evaluation Output**:
   - `network_plot.png`: Gene interaction networks visualization.
   - `volcano_plot.png`: gRNA-gene pairs analysis.
   - IGV files (`.bedgraph` and `bedpe`): Genome browser visualization files.

2. **Analysis Figures**:
   - `knee_plot_scRNA.png`: Knee plot of UMI counts vs. barcode index.
   - `scatterplot_scrna.png`: Scatterplot of total counts vs. genes detected, colored by mitochondrial content.
   - `violin_plot.png`: Distribution of gene counts, total counts, and mitochondrial content.
   - `scRNA_barcodes_UMI_thresholds.png`: Number of scRNA barcodes using different Total UMI thresholds.
   - `guides_per_cell_histogram.png`: Histogram of guides per cell.
   - `cells_per_guide_histogram.png`: Histogram of cells per guide.
   - `guides_UMI_thresholds.png`: Simulating the final number of cells with assigned guides using different minimal number thresholds (at least one guide > threshold value). (Use it to inspect how many cells would have assigned guides. This can be used to check if the final number of cells with guides fit with your expected number of cells)
   - `guides_UMI_thresholds.png`: Histogram of the number of sgRNA represented per cell
   - `cells_per_htp_barplot.png`: Number of Cells across Different HTOs
   - `umap_hto.png`: UMAP Clustering of Cells Based on HTOs (The dimensions represent the distribution of HTOs in each cell)
   - `umap_hto_singlets.png`: UMAP Clustering of Cells Based on HTOs (multiplets removed)

3. **seqSpec Plots**:

   - `seqSpec_check_plots.png`: The frequency of each nucleotides along the Read 1 (Use to inspect the expected read parts with their expected signature) and Read 2 (Use to inspect the expected read parts with their expected signature).

**Structure:**
```
pipeline_dashboard/
  ├── dashboard.html                         
  │
  ├── evaluation_output/                      
  │   ├── network_plot.png                   
  │   ├── volcano_plot.png                  
  │   ├── igv.bedgraph                     
  │   └── igv.bedpe                         
  │
  ├── figures/
  │   ├── knee_plot_scRNA.png                
  │   ├── scatterplot_scrna.png              
  │   ├── violin_plot.png                    
  │   ├── scRNA_barcodes_UMI_thresholds.png  
  │   ├── guides_per_cell_histogram.png      
  │   ├── cells_per_guide_histogram.png      
  │   ├── guides_UMI_thresholds.png          
  │   ├── cells_per_htp_barplot.png          
  │   ├── umap_hto.png                       
  │   └── umap_hto_singlets.png              
  │
  ├── guide_seqSpec_plots/
  │   └── seqSpec_check_plots.png            
  │
  └── hashing_seqSpec_plots/
      └── seqSpec_check_plots.png             
```

### Pipeline metadata

`pipeline_info/` contains run metadata for reproducibility, including the resolved `nextflow.config`, `nextflow.log`, timestamped `params_*.json`, software versions, the original samplesheet copied as `original_samplesheet.csv` or `original_samplesheet.tsv`, and Nextflow execution resource reports. If the samplesheet path comes from a profile or config file, the copied file is taken from that resolved `params.input` value. Demo runs also contain the selected `demo_samplesheet.csv` or `.tsv` and `DEMO_MODE_WARNING.txt`; the warning is duplicated at the output root so reduced-data outputs cannot be mistaken for a full run.

The resource reports are:

| File | Description |
|---|---|
| `execution_trace.txt` | Per-task TSV with process name, status, requested CPUs/memory/time, attempts, wall time, CPU usage, RSS/VMEM, peak RSS/VMEM, I/O counters, container, and work directory. This is the best source for per-process memory and runtime auditing. |
| `execution_report.html` | Nextflow HTML execution report summarizing resource usage across tasks. |
| `execution_timeline.html` | Nextflow HTML task timeline for reviewing task start/end times and concurrency. |

## Pipeline Testing Guide

To ensure proper pipeline functionality, we provide two extensively validated datasets for testing purposes.

### Available Test Datasets

#### 1. TF_Perturb_Seq_Pilot Dataset (Gary-Hon Lab)

The TF_Perturb_Seq_Pilot dataset was generated by the Gary-Hon Lab and is available through the IGVF Data Portal under Analysis Set ID: IGVFDS4389OUWU. To generate the per-sample input file and download the associated FASTQ/configuration files, use the maintained portal download utilities in `download_development/`:

1. First, register for an account on the IGVF Data Portal to obtain your access credentials.

2. Create an IGVF keypair JSON file:

   ```json
   {
     "key": "YOUR_ACCESS_KEY",
     "secret": "YOUR_SECRET_KEY"
   }
   ```

3. Generate the per-sample metadata TSV:

   ```bash
   cd download_development
   pip install -r requirements.txt

   python3 generate_per_sample.py \
       --keypair igvf_key.json \
       --accession IGVFDS4389OUWU \
       --output per_sample.tsv
   ```

4. Download and verify the files, producing a samplesheet with resolved local paths:

   ```bash
   python3 download_igvf.py \
       --sample per_sample.tsv \
       --keypair igvf_key.json \
       --gunzip
   ```

See `download_development/README.md` for optional Google Cloud Storage upload support and fallback seqspec arguments.

All other required input files for running the pipeline with this dataset are already included in the repository under the `example_data` directory.

#### 2. Gasperini et al. Dataset

This dataset comes from a large-scale CRISPR screen study published in Cell ([Gasperini et al., 2019](https://www.cell.com/cell/fulltext/S0092-8674(18)31554-X): "A Genome-wide Framework for Mapping Gene Regulation via Cellular Genetic Screens") and provides an excellent resource for testing the pipeline. The full dataset, including raw sequencing data and processed files, is publicly available through [GEO under accession number GSE120861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861).

### Step-by-Step Testing Instructions

1. **Environment Setup**
   ```bash
   # Clone and enter the repository
   git clone https://github.com/pinellolab/CRISPR_Pipeline.git
   cd CRISPR_Pipeline
   ```

2. **Choose Your Dataset and Follow the Corresponding Instructions:**

   #### Option A: TF_Perturb_Seq_Pilot Dataset
   ```bash
   # Run with LOCAL
   nextflow run main.nf \
      -profile local \
      --input samplesheet.tsv \
      --outdir ./outputs/

   # Run with SLURM
   nextflow run main.nf \
      -profile slurm \
      --input samplesheet.tsv \
      --outdir ./outputs/
   
   # Run with GCP
   nextflow run main.nf \
      -profile google \
      --input samplesheet.tsv \
      --outdir gs://igvf-pertub-seq-pipeline-data/scratch/ # Path to your GCP bucket
   ```

   #### Option B: Gasperini Dataset

   1. Set up the configuration files:
   
      ```bash
      # Copy configuration files and example data
      cp example_gasperini/nextflow.config nextflow.config
      cp -r example_gasperini/example_data/* example_data/
      ```

   2. Obtain sequencing data:
      - Download a subset of the dataset gasperini in your own server.
      - Place files in `example_data/fastq_files` directory

      ```
      NTHREADS=16
      wget https://github.com/10XGenomics/bamtofastq/releases/download/v1.4.1/bamtofastq_linux; chmod +x bamtofastq_linux
      wget https://sra-pub-src-1.s3.amazonaws.com/SRR7967488/pilot_highmoi_screen.1_CGTTACCG.grna.bam.1;mv pilot_highmoi_screen.1_CGTTACCG.grna.bam.1 pilot_highmoi_screen.1_CGTTACCG.grna.bam
      ./bamtofastq_linux --nthreads="$NTHREADS" pilot_highmoi_screen.1_CGTTACCG.grna.bam bam_pilot_guide_1

      wget https://sra-pub-src-1.s3.amazonaws.com/SRR7967482/pilot_highmoi_screen.1_SI_GA_G1.bam.1;mv pilot_highmoi_screen.1_SI_GA_G1.bam.1 pilot_highmoi_screen.1_SI_GA_G1.bam
      ./bamtofastq_linux --nthreads="$NTHREADS" pilot_highmoi_screen.1_SI_GA_G1.bam bam_pilot_scrna_1
      ```
      Now you should see the `bam_pilot_guide_1` and `bam_pilot_scrna_1` directories inside the `example_data/fastq_files` directory. Inside `bam_pilot_guide_1` and `bam_pilot_scrna_1`, there are multiple sets of FASTQ files.

   3. Prepare the whitelist:
      ```bash
      # Extract the compressed whitelist file
      unzip example_data/yaml_files/3M-february-2018.txt.zip
      ```
      Now you should see `3M-february-2018.txt` inside `example_data/yaml_files/` directory.

   4. Launch the pipeline:
      ```bash
      # Run with LOCAL
      nextflow run main.nf \
         -profile local \
         --input samplesheet.tsv \
         --outdir ./outputs/

      # Run with SLURM
      nextflow run main.nf \
         -profile slurm \
         --input samplesheet.tsv \
         --outdir ./outputs/
      
      # Run with GCP
      nextflow run main.nf \
         -profile google \
         --input samplesheet.tsv \
         --outdir gs://igvf-pertub-seq-pipeline-data/scratch/ # Path to your GCP bucket
      ```

### Expected Outputs
The pipeline generates these outputs upon completion:
- `pipeline_outputs`: Contains the final MuData file and cis/trans result tables
- `pipeline_dashboard`: Houses interactive visualization reports and supporting assets only
- `pipeline_dashboard.tar.gz`: Compressed archive of `pipeline_dashboard`

### Troubleshooting
If you encounter any issues during testing:
1. Review log files and intermediate results in the `work/` directory
2. Verify that all input files meet the required format specifications

For additional support or questions, please open an issue on our GitHub repository.

## Credits

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#fg-crispr` channel]

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/crispr for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
