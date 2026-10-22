# CRISPR Pipeline documentation

From raw RNA/guide/HTO reads to auditable perturbation results.

This website documents the current `dev` pipeline: modality-aware input preparation, per-measurement-set QC, guide calling, checkpoint and external-MuData continuation, SCEPTRE/PerTurbo inference, non-targeting calibration, final MuData, and the live execution dashboard.

| If you want to… | Start with… |
|---|---|
| Prepare a new experiment | [Samplesheet](samplesheet.md), [SeqSpec](seqspec.md), [guide metadata](guide-metadata.md) |
| Run end-to-end | [Start here](usage.md) and [run recipes](run-recipes.md) |
| Enable the integrated live dashboard | [Install the bundled plugin](plugin-installation.md) |
| Stop before inference, filter cells, add covariates | [Pause, edit, continue](preinference_checkpoint_user_guide.md) |
| Understand filters or method choices | [QC order](filtering_workflow.md), [covariates](samplesheet_covariates.md), [inference](inference.md) |
| Interpret the outputs | [Output files](output.md), [all MuData fields](mudata_schema.md), [all parameters](parameters.md) |
| See real plots and a completed dashboard | [Results gallery](output-gallery.md) and [Gary Hon dashboard package](https://www.synapse.org/Synapse:syn78349990) |
| Understand the tools behind each stage | [Tools and modules](tool-and-module-map.md) |

The source is [pinellolab/CRISPR_Pipeline](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev). Documentation changes on `dev` are built and deployed by GitHub Actions. The field and parameter references are generated from the checked-in catalog and schema, so changes can be reviewed alongside code. Cite the exact commit recorded in each run's `pipeline_info` when reporting results.

**Plugin note:** a fresh clone does not by itself install the optional live-dashboard Nextflow plugin. Follow the [plugin installation guide](plugin-installation.md) before enabling it. The scientific workflow still runs without online W&B upload.
