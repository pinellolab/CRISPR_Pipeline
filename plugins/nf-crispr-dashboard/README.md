# CRISPR Pipeline live-dashboard plugin

`nf-crispr-dashboard` is a Nextflow trace observer for the detailed CRISPR
Pipeline execution dashboard. It starts the existing self-contained renderer
when a run begins, records live and cached task events, finalizes the HTML
when the run ends, and optionally serves it through a local HTTP port.

On pipeline revisions containing the non-targeting calibration check (PR
#103), its summary, QQ and p-value distribution plots, matched-method
scatterplots, and precision-recall curves appear in the **Inference** block
as the calibration process publishes them. A skipped check shows its reason
there instead. The calibration report itself remains under
`pipeline_outputs/calibration/`.

The plugin is packaged with this pipeline source and is pinned in the pipeline
configuration. Online installations can obtain the plugin from the Nextflow
Plugin Registry once it is published; offline clusters must pre-stage this
exact archive in `NXF_PLUGINS_DIR`. Merely cloning the pipeline source is not
enough for Nextflow to load an uninstalled plugin archive.

For offline installation from the archive bundled in this repository:

```bash
export NXF_PLUGINS_DIR=/shared/path/nextflow-plugins
bash run_support/install_dashboard_plugin.sh
```

Then add `-c conf/dashboard_plugin.config` to the normal `nextflow run`
command. This is an opt-in configuration until the plugin is published to the
registry; enabling an unpublished plugin by default would make new online
installs fail before the analysis begins.

The renderer runs inside the pinned pipeline image by default. Set
`params.DASHBOARD_PLUGIN.rendererCommand` to a list of command tokens when the execution
host uses another container engine or has a pre-pulled image. This command
must end in the container's `python` executable; the plugin appends the
renderer script path and its arguments. No host Python packages are required.
Set `params.DASHBOARD_PLUGIN.imageLibraryDir` to the directory containing a
pre-pulled Nextflow-style `.img` file if that image is not in the configured
Apptainer/Singularity library. Otherwise the plugin uses the pinned OCI URI.

```groovy
params.DASHBOARD_PLUGIN = [
    enabled: true,
    upload: true,                 // optional W&B; falls back to local HTML
    tokenEnv: 'WB_IGVF',
    server: [
        enabled: true,
        host: '127.0.0.1',        // safe default; use SSH port forwarding
        port: 8787,
        refreshSeconds: 30
    ]
]
```

While Nextflow runs, open `http://127.0.0.1:8787/` on the execution host,
or forward the port over SSH. After the run, the same self-contained HTML
remains under `results/pipeline_dashboard/execution_dashboard.html`; the
plugin's in-process HTTP listener stops with Nextflow. W&B upload, if enabled
and configured, also retains the final view online.

Binding the server to a non-loopback address requires the environment
variable named by `params.DASHBOARD_PLUGIN.server.tokenEnv` (default
`CRISPR_DASHBOARD_TOKEN`). HTTP Basic authentication then uses username
`dashboard` and that token as the password. Use TLS termination at a trusted
reverse proxy for any remote exposure; raw HTTP does not encrypt passwords.

Build with `./gradlew assemble` using JDK 21. The distributable ZIP is in
`build/distributions/`. Run `./gradlew installPlugin` for local testing.
