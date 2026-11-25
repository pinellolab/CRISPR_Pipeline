# Optional monitoring and troubleshooting

The supported execution dashboard writes local HTML as the run progresses.
See the [plugin installation guide](plugin-installation.md) for the tested
Nextflow version and the [execution dashboard guide](live_wandb_execution_dashboard.md)
for optional Weights & Biases upload and HTTP serving.

## Axiom telemetry

The Axiom wrapper sends compact workflow and task events. It uses the
`AXIOM_IGVF` environment variable for the token and targets the
`crispr-pipeline` dataset by default. Set the token in your launch environment.

```bash
export AXIOM_OUTDIR=/shared/results/study
export AXIOM_RUN_NAME=study
bin/run_with_axiom.sh nextflow run main.nf -profile slurm \
  --input /shared/inputs/samplesheet.csv --outdir "$AXIOM_OUTDIR" \
  -name "$AXIOM_RUN_NAME"
```

See the [Axiom config](https://github.com/pinellolab/CRISPR_Pipeline/blob/main/conf/axiom.config)
for the upload budget, dashboard refresh and token environment settings.
`AXIOM_INGEST_URL` can select the appropriate regional ingest endpoint.
Dashboard creation requires a token with dashboard permissions. The wrapper
reports upload failures separately from the Nextflow exit status. PNG, SVG
and HTML artifacts remain in the results directory.

## Seqera Platform

Configure the launch environment and credentials in your Seqera Platform
compute environment when tracking there. Keep credentials in the environment
or the platform's secrets store. The local dashboard and pipeline provenance
remain available in the output directory.

## Troubleshoot a run

1. Check the Nextflow launch log for the failed process, attempt and work directory.
2. Open that task's `.command.err`, `.command.out` and `.command.sh`. Check the
   container path and scheduler request against the settings of the compute node.
3. Verify that every node can read the samplesheet paths and container cache.
   For mapping or guide errors, verify the chemistry, read orientation and
   expected modality using [SeqSpec](seqspec.md) and [guide-mapping QC](guide_mapping_overlap_qc.md).
4. Check `pipeline_info/` for the source commit, resolved config, inputs and
   trace. Use the same source, launch directory and work directory for
   [task resume](run-recipes.md#retry-a-failed-invocation).
