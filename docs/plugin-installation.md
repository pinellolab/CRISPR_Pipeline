# Install the live-dashboard plugin

The integrated dashboard uses the `nf-crispr-dashboard@0.1.0` Nextflow plugin.
For now, the versioned plugin archive is **bundled in this repository**, but
Nextflow does not load it directly from the source tree. Install it once into
an absolute, writable `NXF_PLUGINS_DIR` visible to the Nextflow launcher before
adding `conf/dashboard_plugin.config` to a run. This is optional: the analysis
can run without the plugin, but a direct Nextflow invocation then does not
provide this integrated live dashboard.

## One-time installation on each execution environment

From the root of the same pinned pipeline checkout that will run the analysis:

```bash
cd /absolute/path/CRISPR_Pipeline
export NXF_PLUGINS_DIR=/absolute/shared/path/nextflow-plugins
bash run_support/install_dashboard_plugin.sh
```

Use a persistent directory with enough space and write permission. On a
cluster, make it accessible from the node that launches Nextflow; use the same
`NXF_PLUGINS_DIR` for later resumes. On the Pinello local server, choose the
approved scratch plugin cache, **not** your home directory or system `/tmp`.
The installer verifies the bundled ZIP's pinned SHA-256, refuses to overwrite
a different installed version, and is safe to rerun for the same archive.

Verify the installation without starting a pipeline:

```bash
test -f "$NXF_PLUGINS_DIR/nf-crispr-dashboard-0.1.0/.archive.sha256"
sha256sum plugins/nf-crispr-dashboard/dist/nf-crispr-dashboard-0.1.0.zip
```

The archive SHA-256 must be
`6fb2b20f5012d65b2ed1e1e14c026d49f1f4a779faee1d20a2fbac870b9919b0`.
The installer prints `Installed dashboard plugin` or `Dashboard plugin already
installed` on success. The actual version and checksum are defined in
`run_support/install_dashboard_plugin.sh`; if they change, use the matching
checkout and reinstall rather than copying an older plugin directory.

## Enable it for a run

Keep `NXF_PLUGINS_DIR` exported in the shell that starts Nextflow. Add the
plugin config after your dataset/resource configs and before the generated
provenance config, for example:

```bash
nextflow \
  -c /absolute/path/dataset.resources.config \
  -c /absolute/path/CRISPR_Pipeline/conf/dashboard_plugin.config \
  run /absolute/path/CRISPR_Pipeline/main.nf \
  -params-file /absolute/path/dataset.params.json \
  --input /absolute/path/samplesheet.csv \
  --outdir /absolute/path/results \
  --ANALYSIS_mode full \
  -work-dir /absolute/path/work \
  -name example_full_01
```

The example omits site-specific executor, container-storage, and provenance
configs; follow the [run recipes](run-recipes.md) for those. On the Pinello
local server, source `run_support/container_storage.sh` **immediately before
every launch or resume**, check its return status, and add
`run_support/container_storage.config` after dataset resources and before
`provenance.generated.config`. Record the plugin config, installer, archive,
storage files, and run wrapper as provenance artifacts. A remote SLURM or
cloud site needs its own accessible storage paths; do not copy Pinello paths.

The opt-in config enables local HTML generation and attempts W&B upload with
the token environment variable `WB_IGVF`. Do not put a token in a config or
command line. Without a working W&B token, keep or set `upload: false` in a
run-specific config override; local HTML still works. The final self-contained
file is `results/pipeline_dashboard/execution_dashboard.html`. See the
[live dashboard guide](live_wandb_execution_dashboard.md) for optional local
HTTP serving, port forwarding, authentication, and update behavior.

The plugin is **not yet assumed available from the Nextflow Plugin Registry**.
Do not remove this installation step or enable the plugin by default until a
clean machine can resolve the published version. The source and publishing
instructions are in `plugins/nf-crispr-dashboard/README.md`.

## Clean-install verification

On 22 October 2026, a fresh `dev` clone at commit `16a07e4` was tested with
an isolated Nextflow home and plugin directory. The official launcher
installed Nextflow 26.04.6; the bundled installer verified the archive and
installed `nf-crispr-dashboard@0.1.0`; `nextflow config` resolved the plugin
declaration and dashboard settings. No analysis workflow was launched. This
checks reproducible installation and config parsing, **not** runtime observer
startup, container access, or scientific results.
