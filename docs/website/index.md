# CRISPR Pipeline

The CRISPR Pipeline turns RNA, guide and optional hashtag reads from
single-cell perturbation screens into quality-controlled count matrices and
perturbation results. It records what changed at each stage, from SeqSpec
validation and mapping through guide assignment and inference.

The current documented release is **v1.1.0**. Pin that tag for reproducible
runs and report the source commit recorded in `pipeline_info/` when sharing
results. [Release source](https://github.com/pinellolab/CRISPR_Pipeline/tree/v1.1.0)

## Start running

- [Install and prepare inputs](usage.md), including the
  [samplesheet](samplesheet.md), [SeqSpec](seqspec.md), and [guide/HTO metadata](guide-metadata.md).
- [Run on Slurm](run-recipes.md#slurm-configuration) with Nextflow 26.04.6 and
  the [bundled execution dashboard plugin](plugin-installation.md).
- [Download current configs and study examples](examples.md), including
  CC-Perturb-seq, TAP-seq, Gasperini, Replogle, and the Gary Hon TF benchmark.
- [Stop at Post Concatenation QC and continue later](preinference_checkpoint_user_guide.md).
  The default checkpoint publishes the QC report and a MuData handoff.

## Understand results

- [QC and filtering](filtering_workflow.md) and
  [Post Concatenation QC](postconcat_embedding_qc.md) explain when cells enter and leave each stage.
- [Inference and calibration](inference.md) cover local SCEPTRE/PerTurbo tests,
  transcriptome-wide PerTurbo tests, and the non-targeting-control diagnostic.
- [Outputs](output.md), [MuData fields](mudata_schema.md), and the
  [parameter reference](parameters.md) describe what the run publishes.
- [Results gallery](output-gallery.md) shows real figures and dashboards.
  [Temporary benchmark review](temporary-benchmarks.md)

## Project and support

[Source and changelog](https://github.com/pinellolab/CRISPR_Pipeline) ·
[Report an issue](https://github.com/pinellolab/CRISPR_Pipeline/issues) ·
[Tool and method citations](https://github.com/pinellolab/CRISPR_Pipeline/blob/main/CITATIONS.md)

The site is built from the pipeline's `main` branch. Parameter and MuData
references are generated from the checked-in schema and field catalog.
