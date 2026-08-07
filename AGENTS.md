# Agent instructions

## Nextflow run provenance

For every prepared or resumed dataset run:

1. Run `bin/run_provenance_agent.py prepare` immediately before Nextflow.
2. Pass the exact pipeline repository, actual executing source directory, output
   `pipeline_info` directory, run name, and every dataset-specific config,
   params, samplesheet, guide metadata, and wrapper as named `--artifact`
   arguments.
3. Use `--expected-commit` when a run is intentionally pinned to an older source
   snapshot. Otherwise the agent uses the current `HEAD`.
4. Do not use `--allow-dirty` or `--allow-source-mismatch` unless the user
   explicitly accepts non-commit source provenance.
5. Run the agent's `check` command after preparation and before launching
   Nextflow. Stop if the source or any recorded artifact changed.
6. Pass the generated `provenance.generated.config` to Nextflow after other
   `-c` files so the resolved configuration reports the exact commit.
7. Keep `repository_provenance.json`, `repository_provenance.tsv`,
   `pipeline_source_files.sha256`, and `provenance.generated.config` in the
   run's output `pipeline_info` directory.

Never hardcode the repository's current commit in a tracked config file. The
commit changes when that file is committed; generate commit-specific metadata
at run time instead.
