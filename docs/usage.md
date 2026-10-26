# Start here

CRISPR Pipeline maps RNA, guide, and optional HTO reads; performs quality control and guide assignment; then runs local SCEPTRE and local/global PerTurbo inference. Input is a **modality-aware samplesheet**, not the older nf-core `sample,fastq_1,fastq_2` format. Read the [samplesheet contract](samplesheet.md) before launching.

## Prerequisites

- Nextflow and an Apptainer/Singularity or Docker runtime appropriate for your machine.
- A clean, pinned pipeline checkout and persistent work directory. Use absolute input paths; for SLURM, every compute node must see the same work and input paths.
- A complete samplesheet, RNA and guide SeqSpec YAMLs, guide metadata, barcode onlist, and a compatible GTF/genome build. HTO assays also need an HTO SeqSpec and metadata.
- Enough space for references, container images, mapping, and inference. On the Pinello server, follow the [local storage policy](https://github.com/pinellolab/CRISPR_Pipeline/blob/dev/AGENTS.md) and the preparation preflight; other hosts need their own storage paths.

The example input files in this repository explain the format but have historical cloud URLs; replace those with paths accessible to your executor. Check guide sequences, read orientation, barcode chemistry, and coordinate build against **your** assay.

## Choose the run mode

`ANALYSIS_mode` defaults to `checkpoint` in the current `dev` configuration. It completes mapping, QC, and guide assignment, then publishes a pre-inference MuData without running inference. Use `full` for an unattended end-to-end analysis, or `continue` to validate an edited checkpoint and run inference. These are distinct Nextflow invocations. `-resume` only retries cached tasks in an invocation; it does not convert a checkpoint into a continuation.

```bash
nextflow run /absolute/path/CRISPR_Pipeline/main.nf \
  -profile local \
  -c /absolute/path/resources.config \
  --input /absolute/path/samplesheet.csv \
  --outdir /absolute/path/results_full \
  --ANALYSIS_mode full \
  --REFERENCE_gtf_local_path /absolute/path/reference.gtf.gz \
  -work-dir /absolute/path/work_full \
  -name example_full_01
```

The `-profile` value and resource config are site-specific. The `local` profile is suitable only when running directly on one machine; the `slurm` profile requires a site-specific executor/partition config and shared paths. Singularity containers are enabled in the current pipeline config. For a dataset-specific reproducible run, use a params file and the [run recipes](run-recipes.md), including the required provenance and container-storage checks on the Pinello server.

## The first checks

Before paying for a full run, verify that the SeqSpecs extract the intended cell barcode, UMI, cDNA and guide sequence; that guide metadata has unique IDs, exact captured spacer sequences and the correct coordinate build; and that RNA and guide mapping recover overlapping barcodes. `DEMO_MODE=true` selects one complete measurement set for a small pre-run check, but its outputs are **not** full-dataset results.

After launch, inspect the live [execution dashboard](live_wandb_execution_dashboard.md). Each measurement set has its own mapping and preprocessing QC. The final [MuData dictionary](mudata_schema.md), [output guide](output.md) and [parameter reference](parameters.md) explain the resulting data.

For an integrated live dashboard in a direct Nextflow run, first follow the
[bundled plugin installation guide](plugin-installation.md). The plugin is
currently opt-in and is not installed merely by cloning this repository.

## What happens next

1. SeqSpec validation and reference preparation.
2. RNA, guide, and optional HTO mapping, followed by per-measurement-set RNA cell calling and QC.
3. RNA/guide barcode intersection, guide assignment, optional HTO singlet and clone filters, then post-concatenation mitochondrial/gene QC.
4. Inference using the **raw-count** MuData; temporary normalized PCA/UMAP data are visualization-only.
5. Calibration, evaluation, final tables, the dashboard, and provenance artifacts.

See [QC and filtering order](filtering_workflow.md) and [inference](inference.md) for the exact behavior and configurable exceptions.
