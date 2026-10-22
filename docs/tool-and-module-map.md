# Tools and pipeline modules

The pipeline orchestrates tools inside versioned containers; users do not
install their Python/R dependencies into the Nextflow launcher environment.
The container image used by a particular run is recorded in its resolved
config and task trace. These links are to the **upstream projects' official
documentation**; the pipeline's parameters and wrappers determine which
features are actually used here.

| Stage | This repository | Upstream tool and documentation | What the pipeline does |
|---|---|---|---|
| Orchestration and provenance | [`main.nf`](https://github.com/pinellolab/CRISPR_Pipeline/blob/dev/main.nf), [`workflows/crispr_pipeline/main.nf`](https://github.com/pinellolab/CRISPR_Pipeline/blob/dev/workflows/crispr_pipeline/main.nf) | [Nextflow documentation](https://docs.seqera.io/nextflow/) | Splits the samplesheet by modality and measurement set, schedules containerized tasks, caches/resumes, records reports. |
| Read structure | [`seqSpecCheck_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/seqSpecCheck_pipeline) | [seqspec](https://github.com/pachterlab/seqspec) | Checks and parses RNA, guide, and optional HTO library structures before mapping. |
| RNA and guide mapping | [`mapping_rna_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/mapping_rna_pipeline), [`mapping_guide_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/mapping_guide_pipeline) | [kallisto, bustools, and kb-python](https://kallisto.readthedocs.io/en/latest/) | Builds/uses references and counts molecules; the guide branch has base-editing mapping variants. |
| Cell QC and embedding | [`preprocessing_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/preprocessing_pipeline), [`postconcat_embedding_qc`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/modules/local/postconcat_embedding_qc) | [Scanpy preprocessing](https://scanpy.readthedocs.io/en/stable/api/preprocessing.html), [Scrublet in Scanpy](https://scanpy.scverse.org/en/stable/api/generated/scanpy.pp.scrublet.html) | Applies per-set cell QC, then post-concatenation QC/visualization. Normalized PCA/UMAP arrays are temporary; inference receives raw counts. |
| MuData and HTO | [`CreateMuData`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/modules/local/CreateMuData), [`demultiplex`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/modules/local/demultiplex) | [MuData](https://mudata.scverse.org/en/latest/), [GMM-Demux](https://github.com/CHPGenetics/GMM-demux) | Joins modalities by qualified barcodes and optionally demultiplexes HTOs; HTO singlet filtering is conditional on assay configuration. |
| Guide assignment | [`guide_assignment_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/guide_assignment_pipeline) | [SCEPTRE](https://katsevich-lab.github.io/sceptre/), [CLEANSER](https://github.com/Gersbachlab-Bioinformatics/CLEANSER) | Chooses one configured assignment method, then applies guide-positive and per-cell guide-count QC. CLEANSER is optional, not an additional inference method. |
| Perturbation inference | [`inference_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/inference_pipeline) | [SCEPTRE](https://katsevich-lab.github.io/sceptre/), [PerTurbo](https://github.com/pinellolab/PerTurbo) | Local/cis SCEPTRE and PerTurbo, plus global PerTurbo when selected; prepares pairs, control regime, covariates, and merged result tables. |
| Calibration and dashboard | [`ntc_calibration_pipeline`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/subworkflows/local/ntc_calibration_pipeline), [`nf-crispr-dashboard`](https://github.com/pinellolab/CRISPR_Pipeline/tree/dev/plugins/nf-crispr-dashboard) | [W&B documentation](https://docs.wandb.ai/) (optional upload) | Matched non-targeting evaluation when possible; local live/final dashboard from the bundled Nextflow plugin, with optional W&B publication. |

```text
samplesheet + SeqSpecs + metadata
    → RNA / guide / optional HTO mapping (per measurement set)
    → cell QC and MuData assembly → guide assignment / optional HTO QC
    → post-concatenation QC → checkpoint or inference
    → local/global tables + calibration + dashboard + provenance
```

The [filtering guide](postconcat_embedding_qc.md) describes the exact QC order;
the [inference guide](inference.md) describes pair selection and covariates.
Use the [run's provenance](output.md#pipeline-information), not a floating
upstream documentation page, to establish its actual source commit,
containers, and parameters. A method listed here may be skipped by an assay
or run mode.
