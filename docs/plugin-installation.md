# Install the live-dashboard plugin

The integrated dashboard uses the `nf-crispr-dashboard@0.1.2` Nextflow plugin.
For now, the versioned plugin archive is **bundled in this repository**, but
Nextflow does not load it directly from the source tree. Install it once into
an absolute, writable `NXF_PLUGINS_DIR` visible to the Nextflow launcher.
The v1.1.0 root config already includes `conf/dashboard_plugin.config`.
Install the plugin before loading that config. A custom configuration that
omits the observer can run analysis without it when IGVF packaging is disabled.

## One-time installation on each execution environment

From the root of the same pinned pipeline checkout that will run the analysis:

```bash
cd /absolute/path/CRISPR_Pipeline
export NXF_VER=26.04.6
export NXF_SYNTAX_PARSER=v1
export NXF_PLUGINS_DIR=/absolute/shared/path/nextflow-plugins
bash run_support/install_dashboard_plugin.sh
nextflow -version
```

Nextflow 26.04.6 is the engine version used to build and validate the bundled
plugin. This version pin applies to the Nextflow launcher on Slurm and other
executors. Export these variables in the shell or batch script that invokes
Nextflow.

Use a persistent directory with enough space and write permission. On a
cluster, make it accessible from the node that launches Nextflow; use the same
`NXF_PLUGINS_DIR` for later resumes. On the Pinello local server, choose the
approved scratch plugin cache, **not** your home directory or system `/tmp`.
The installer verifies the bundled ZIP's pinned SHA-256, refuses to overwrite
a different installed version, and is safe to rerun for the same archive.

Verify the installation without starting a pipeline:

```bash
test -f "$NXF_PLUGINS_DIR/nf-crispr-dashboard-0.1.2/.archive.sha256"
sha256sum plugins/nf-crispr-dashboard/dist/nf-crispr-dashboard-0.1.2.zip
```

On macOS, use `shasum -a 256` instead of `sha256sum`. The archive SHA-256 must be
`02a1eba42b8f30dc29952a603a36c485a089dfa5ea18d7104e415040ec902bbf`.
The installer prints `Installed dashboard plugin` or `Dashboard plugin already
installed` on success. The actual version and checksum are defined in
`run_support/install_dashboard_plugin.sh`; if they change, use the matching
checkout and reinstall rather than copying an older plugin directory.

The bundled installer uses the prebuilt plugin. For development, see the
[plugin source guide](https://github.com/pinellolab/CRISPR_Pipeline/blob/v1.1.0/plugins/nf-crispr-dashboard/README.md)
for `./gradlew installPlugin`.

## Enable it for a run

Version 0.1.2 records task work directories and attempts, retains observer
history across resumed runs, and supports a pinned local renderer image.
The renderer publishes bounded task diagnostics alongside the execution
history. See the [execution dashboard](live_wandb_execution_dashboard.md) for
archive limits and recovery-status interpretation.

When using Nextflow 26.04.6 with this pipeline's legacy DSL2 closures, export
`NXF_SYNTAX_PARSER=v1` before configuration checks and execution. A preview
can verify graph construction without executing tasks; disable
`OUTPUT_IGVF_PACKAGE` for that preview only, because preview does not produce
the final files required by the completion packaging hook. Full runs retain
the default IGVF packaging and TSV outputs. Completion packaging resolves
workflow metadata through the module method binding, so deferred DSL2 hook
delegation cannot hide the workflow's final status. Keep the launcher log at
`<launchDir>/.nextflow.log` so the completion step can publish the required
`pipeline_info/nextflow.log` into the IGVF package.

Keep `NXF_PLUGINS_DIR` exported in the shell that starts Nextflow. With a
custom config stack, load the plugin config before the generated
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

The release config enables local HTML generation and attempts W&B upload with
the token environment variable `WB_IGVF`. Do not put a token in a config or
command line. Without a working W&B token, keep or set `upload: false` in a
run-specific config override; local HTML still works. The final self-contained
file is `results/pipeline_dashboard/execution_dashboard.html`. See the
[live dashboard guide](live_wandb_execution_dashboard.md) for optional local
HTTP serving, port forwarding, authentication, and update behavior.

The release uses the bundled plugin archive for installation. The source and publishing
instructions are in `plugins/nf-crispr-dashboard/README.md`.

## Slurm renderer troubleshooting

The observer starts its renderer on the Nextflow launcher host. Its default
renderer command uses Apptainer. If `pipeline_dashboard/.live/` exists but
`execution_dashboard.html` is missing, read
`pipeline_dashboard/.live/monitor.log` and verify the launcher can invoke the
renderer and access the output and work directories.

For a local image, use an absolute path such as
`params.DASHBOARD_PLUGIN.image = '/shared/images/crispr-base.sif'`. The
observer handles absolute image paths directly. A site with only Singularity
can set `params.DASHBOARD_PLUGIN.rendererCommand` to a command prefix that
runs Python in the compatible image and binds the required paths. See the
[dashboard guide](live_wandb_execution_dashboard.md) for renderer settings.

## Clean-install verification

On 22 October 2026, a fresh `dev` clone at commit `16a07e4` was tested with
an isolated Nextflow home and plugin directory. The official launcher
installed Nextflow 26.04.6; the bundled installer verified the archive and
installed `nf-crispr-dashboard@0.1.0`; `nextflow config` resolved the plugin
declaration and dashboard settings. No analysis workflow was launched. This
checks reproducible installation and config parsing, **not** runtime observer
startup, container access, or scientific results.
