<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-crispr_logo_dark.png">
    <img alt="nf-core/crispr" src="docs/images/nf-core-crispr_logo_light.png">
  </picture>
</h1>

[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/crispr)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23crispr-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/crispr)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

# CRISPR Pipeline

A comprehensive pipeline for single-cell Perturb-Seq analysis that enables robust processing and analysis of CRISPR screening data at single-cell resolution.

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
- `pairs_to_test.csv`: Defines perturbation pairs for comparison analysis (required if testing predefined pairs)

For detailed specifications, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

## Running the Pipeline 

### Pipeline Configuration

Before running the pipeline, customize the configuration files for your environment:

#### 1. Data and Analysis Parameters (`nextflow.config`)

Update the pipeline-specific parameters in the `params` section, for example:

```groovy
// Input data paths
input = "/path/to/your/samplesheet.csv"

// Reference data (update paths as needed)
METADATA_sgRNA = "/path/to/guide_metadata.tsv" 
METADATA_hash = "/path/to/hash_metadata.tsv"
SEQUENCE_PARSING_barcode_list = "/path/to/barcode_list.txt"

// Analysis parameters (adjust for your experiment)
QC_min_genes_per_cell = 500
QC_min_cells_per_gene = 3
QC_pct_mito = 20
GUIDE_ASSIGNMENT_method = 'sceptre'  // or 'cleanser'
INFERENCE_method = 'perturbo'        // or 'sceptre'
```

#### 2. Compute Environment Configuration

Choose and configure your compute profile by updating the relevant sections:

##### 🖥️ **Local Development**
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

// Run with: nextflow run main.nf -profile google
```

#### 3. Container Configuration

The pipeline uses pre-built containers. Update if you have custom versions:

```groovy
containers {
    base     = 'sjiang9/conda-docker:0.2'      // Main analysis tools
    cleanser = 'sjiang9/cleanser:0.3'          // Guide assignment
    sceptre  = 'sjiang9/sceptre-igvf:0.1'     // Guide assignment/Perturbation inference 
    perturbo = 'loganblaine/perturbo:latest'   // Perturbation inference 
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

The output files will be generated in the `pipeline_outputs` and `pipeline_dashboard` directory.

### Generated Files

Within the `pipeline_outputs` directory, you will find:

- inference_mudata.h5mu - MuData format output
- per_element_output.tsv - Per-element analysis
- per_guide_output.tsv - Per-guide analysis

**Structure:**

```
📁 pipeline_outputs/
   ├── 📄 inference_mudata.h5mu    
   ├── 📄 per_element_output.tsv    
   └── 📄 per_guide_output.tsv     
```

For details, see our [documentation](https://docs.google.com/document/d/1Z1SOlekIE5uGyXW41XxnszxaYdSw0wdAOUVzfy3fj3M/edit?tab=t.0#heading=h.ctbx1w9hj619).

### Generated Figures

The pipeline produces several figures:

Within the `pipeline_dashboard` directory, you will find:

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
📁 pipeline_dashboard/
  ├── 📄 dashboard.html                         
  │
  ├── 📁 evaluation_output/                      
  │   ├── 🖼️ network_plot.png                   
  │   ├── 🖼️ volcano_plot.png                  
  │   ├── 📄 igv.bedgraph                     
  │   └── 📄 igv.bedpe                         
  │
  ├── 📁 figures/
  │   ├── 🖼️ knee_plot_scRNA.png                
  │   ├── 🖼️ scatterplot_scrna.png              
  │   ├── 🖼️ violin_plot.png                    
  │   ├── 🖼️ scRNA_barcodes_UMI_thresholds.png  
  │   ├── 🖼️ guides_per_cell_histogram.png      
  │   ├── 🖼️ cells_per_guide_histogram.png      
  │   ├── 🖼️ guides_UMI_thresholds.png          
  │   ├── 🖼️ cells_per_htp_barplot.png          
  │   ├── 🖼️ umap_hto.png                       
  │   └── 🖼️ umap_hto_singlets.png              
  │
  ├── 📁 guide_seqSpec_plots/
  │   └── 🖼️ seqSpec_check_plots.png            
  │
  └── 📁 hashing_seqSpec_plots/
      └── 🖼️ seqSpec_check_plots.png             
```

## Pipeline Testing Guide

To ensure proper pipeline functionality, we provide two extensively validated datasets for testing purposes.

### Available Test Datasets

#### 1. TF_Perturb_Seq_Pilot Dataset (Gary-Hon Lab)

The TF_Perturb_Seq_Pilot dataset was generated by the Gary-Hon Lab and is available through the IGVF Data Portal under Analysis Set ID: IGVFDS4389OUWU. To access the fastq files, you need to:

1. First, register for an account on the IGVF Data Portal to obtain your access credentials.

2. Once you have your credentials, you can use our provided Python script to download all necessary FASTQ files:

   ```bash
   cd example_data
   python download_fastq.py \
       --sample per-sample_file.tsv \
       --access-key YOUR_ACCESS_KEY \
       --secret-key YOUR_SECRET_KEY
   ```
   
   💡 **Note:** You'll need to replace `YOUR_ACCESS_KEY` and `YOUR_SECRET_KEY` with the credentials from your IGVF portal account. These credentials can be found in your IGVF portal profile settings.

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
The pipeline generates two directories upon completion:
- `pipeline_outputs`: Contains all analysis results
- `pipeline_dashboard`: Houses interactive visualization reports

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

