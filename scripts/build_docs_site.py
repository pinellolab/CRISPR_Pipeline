#!/usr/bin/env python3
"""Assemble the GitHub Pages source from canonical checked-in documentation.

No pipeline environment or network connection is needed for this step.
"""

import csv
import json
import re
import shutil
from html import escape
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from render_slurm_config import CONTAINER_BLOCK, PARAM_BLOCK, assignments


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DEST = ROOT / "build" / "docs-site"

PAGES = (
    "temporary-benchmarks.md",
    "usage.md",
    "plugin-installation.md",
    "configuration.md",
    "samplesheet.md",
    "seqspec.md",
    "guide-metadata.md",
    "samplesheet_covariates.md",
    "run-recipes.md",
    "telemetry.md",
    "preinference_checkpoint_user_guide.md",
    "filtering_workflow.md",
    "postconcat_embedding_qc.md",
    "inference.md",
    "live_wandb_execution_dashboard.md",
    "tool-and-module-map.md",
    "output.md",
    "output-gallery.md",
    "mudata_schema.md",
    "qc_outputs.md",
    "target_aware_gene_filter.md",
    "guide_mapping_overlap_qc.md",
    "clone_and_saturation_qc.md",
    "examples.md",
)


DATASET_EXAMPLES = (
    "gasperini_at_scale",
    "replogle_k562_essential",
    "gary_hon_tf_benchmark",
)

ENGINE_DESCRIPTIONS = {
    "max_cpus": "Upper CPU limit used by the pipeline's task resource selectors.",
    "max_memory": "Upper memory limit used by the pipeline's task resource selectors.",
    "INFERENCE_SCEPTRE_formula_object": "Legacy SCEPTRE setting. Use the shared inference covariate options for the current model.",
    "INFERENCE_SCEPTRE_WRITE_MUDATA": "Write intermediate MuData containing SCEPTRE results.",
    "INFERENCE_PERTURBO_RESULT_BATCH_ROWS": "Rows per batch during global PerTurbo result conversion.",
    "INFERENCE_PERTURBO_MAX_BH_WORKING_BYTES": "Host memory budget in bytes for Benjamini–Hochberg correction.",
    "INFERENCE_PERTURBO_COMPACT_RESULT_FLOATS": "Use compact floating-point columns when converting PerTurbo results.",
    "publish_dir_mode": "Nextflow publication mode for output files.",
    "css": "Dashboard stylesheet asset directory.",
    "js": "Dashboard JavaScript asset directory.",
    "svg": "Dashboard SVG asset directory.",
    "DEBUG_VAR": "Enable diagnostic logging.",
    "google_bucket": "Google Cloud Storage bucket used by the Google Batch profile.",
    "google_project": "Google Cloud project for Google Batch execution.",
    "google_region": "Region for Google Batch execution.",
    "email": "Email address for run completion notifications.",
    "email_on_fail": "Email address for failure notifications.",
    "plaintext_email": "Send run notification emails in plain text.",
    "monochrome_logs": "Disable colored pipeline log messages.",
    "hook_url": "Optional webhook URL for run notifications.",
    "validate_params": "Enable pipeline parameter validation.",
    "version": "Print the pipeline version.",
    "pipelines_testdata_base_path": "Base URL for pipeline test data.",
    "config_profile_name": "Profile name shown in run metadata.",
    "config_profile_description": "Profile description shown in run metadata.",
    "custom_config_version": "Revision used to construct the nf-core config URL.",
    "custom_config_base": "Base URL for nf-core config files.",
    "config_profile_contact": "Contact shown in profile metadata.",
    "config_profile_url": "Website shown in profile metadata.",
    "DASHBOARD_PLUGIN.enabled": "Enable the execution dashboard observer.",
    "DASHBOARD_PLUGIN.upload": "Enable optional Weights & Biases upload.",
    "DASHBOARD_PLUGIN.project": "Weights & Biases project for uploaded runs.",
    "DASHBOARD_PLUGIN.tokenEnv": "Name of the environment variable for the upload token.",
    "DASHBOARD_PLUGIN.pollSeconds": "Observer polling interval in seconds.",
    "DASHBOARD_PLUGIN.server.enabled": "Serve the live dashboard over HTTP.",
    "DASHBOARD_PLUGIN.server.host": "Dashboard server bind address.",
    "DASHBOARD_PLUGIN.server.port": "Dashboard server port.",
    "DASHBOARD_PLUGIN.server.refreshSeconds": "Dashboard refresh interval in seconds.",
    "DASHBOARD_PLUGIN.server.tokenEnv": "Name of the environment variable for the dashboard server token.",
}


def publish_dataset_examples(destination):
    for name in DATASET_EXAMPLES:
        source = ROOT / "examples/datasets" / name
        target = destination / name
        shutil.copytree(source, target, dirs_exist_ok=True)
        shutil.copyfile(ROOT / "examples/prepare_samplesheet.py", target / "prepare_samplesheet.py")
        with ZipFile(destination / f"{name}.zip", "w", compression=ZIP_DEFLATED, strict_timestamps=False) as archive:
            for path in sorted(target.rglob("*")):
                if path.is_file():
                    archive.write(path, path.relative_to(destination))


def flatten_properties(node, prefix=""):
    for key, value in node.get("properties", {}).items():
        name = f"{prefix}.{key}" if prefix else key
        if "properties" in value:
            yield from flatten_properties(value, name)
        else:
            yield name, value


def render_parameter_table(rows, specs):
    lines = [
        "| Parameter | Default on `main` | Description |",
        "|---|---|---|",
    ]
    for name, value in rows:
        spec = specs.get(name, {})
        description = spec.get("description") or ENGINE_DESCRIPTIONS.get(name)
        if name.startswith("containers."):
            description = f"Container image for {name.removeprefix('containers.')} tasks."
        if not description:
            raise ValueError(f"Missing configuration description for {name}")
        if options := spec.get("enum"):
            description += " Allowed values are " + ", ".join(f"`{option}`" for option in options) + "."
        if spec:
            types = spec.get("type", "unspecified")
            type_text = " or ".join(types) if isinstance(types, list) else types
            details = f"Schema type is {type_text}. {spec.get('help_text', '')}".strip()
            description += f" <details><summary>Details</summary>{escape(details)}</details>"
        description = " ".join(description.split()).replace("|", "&#124;")
        default = f"`{value.replace('|', '&#124;')}`" if value is not None else "—"
        lines.append(f"| `{name}` | {default} | {description} |")
    return "\n".join(lines)


def parameter_groups(block):
    for paragraph in re.split(r"\n\s*\n", block.strip("\n")):
        rows = list(assignments(paragraph))
        if rows:
            heading = paragraph.splitlines()[0].strip().removeprefix("// ")
            yield heading, rows


def dashboard_defaults(config):
    block = config.split("params.DASHBOARD_PLUGIN = [\n", 1)[1].split("\n]", 1)[0]
    path = ["DASHBOARD_PLUGIN"]
    for line in block.splitlines():
        line = line.strip().removesuffix(",")
        if line == "]":
            path.pop()
        elif ": " in line and not line.startswith("//"):
            name, value = line.split(": ", 1)
            if value == "[":
                path.append(name)
            else:
                yield ".".join([*path, name]), value


def render_configuration():
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    specs = {name: spec for section in schema["$defs"].values() for name, spec in flatten_properties(section)}
    basic, *advanced = PARAM_BLOCK.findall((ROOT / "nextflow.config").read_text())
    advanced = "\n\n".join(advanced)
    containers = CONTAINER_BLOCK.search(advanced).group(1)
    dashboard = (ROOT / "conf/dashboard_plugin.config").read_text()
    basic_groups = list(parameter_groups(basic))
    advanced_groups = [
        *parameter_groups(CONTAINER_BLOCK.sub("", advanced)),
        ("Runtime images", [(f"containers.{key}", value) for key, value in assignments(containers)]),
        ("Dashboard observer", list(dashboard_defaults(dashboard))),
    ]
    configured = {name for _, rows in basic_groups + advanced_groups for name, _ in rows}
    if schema_only := [name for name in specs if name not in configured]:
        advanced_groups.append(("Schema-only parameters", [(name, None) for name in schema_only]))

    page = (DOCS / "configuration.md").read_text()
    for label, groups in (("BASIC", basic_groups), ("ADVANCED", advanced_groups)):
        marker = f"<!-- {label}_PARAMETERS -->"
        if page.count(marker) != 1:
            raise ValueError(f"Expected one {marker} in configuration.md")
        tables = "\n\n".join(f"### {title}\n\n{render_parameter_table(rows, specs)}" for title, rows in groups)
        page = page.replace(marker, tables)
    return page


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for name in PAGES:
        if name != "configuration.md":
            shutil.copyfile(DOCS / name, DEST / name)
    (DEST / "configuration.md").write_text(render_configuration())
    shutil.copyfile(DOCS / "website" / "index.md", DEST / "index.md")
    styles = DEST / "stylesheets"
    styles.mkdir(exist_ok=True)
    shutil.copyfile(DOCS / "website" / "site.css", styles / "site.css")
    shutil.copytree(DOCS / "images" / "result_examples", DEST / "images" / "result_examples", dirs_exist_ok=True)
    shutil.copyfile(DOCS / "mudata_schema_catalog.tsv", DEST / "mudata_schema_catalog.tsv")
    (DEST / "parameters.md").write_text(
        "# Parameter reference\n\nSee [Basic and advanced configuration](configuration.md) for every parameter.\n"
    )

    examples = DEST / "downloadable-examples"
    examples.mkdir(exist_ok=True)
    for source, target in (
        ("nextflow_cc.config", "nextflow_cc.config"),
        ("nextflow_tapseq.config", "nextflow_tapseq.config"),
        ("conf/simplified_pipeline/nextflow_yang.dev.config", "nextflow_yang.dev.config"),
        ("conf/simplified_pipeline/nextflow_slurm_all_options.config", "nextflow_slurm_all_options.config"),
        ("docs/examples/minimal_samplesheet.csv", "samplesheet.csv"),
        ("docs/examples/hashing_samplesheet.csv", "hashing_samplesheet.csv"),
        ("example_data/guide_metadata.tsv", "guide_metadata.tsv"),
        ("example_data/hash_metadata.tsv", "hash_metadata.tsv"),
        ("example_data/yaml_files/rna_seqspec.yml", "rna_seqspec.yml"),
        ("example_data/yaml_files/guide_seqspec.yml", "guide_seqspec.yml"),
        ("example_data/yaml_files/hash_seqspec.yml", "hash_seqspec.yml"),
    ):
        shutil.copyfile(ROOT / source, examples / target)
    publish_dataset_examples(examples)
    # Validate that the authoritative field catalog is parseable and has no
    # duplicate path/field entries before publishing its rendered page.
    with (DOCS / "mudata_schema_catalog.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [(row["path"], row["field"]) for row in rows]
    if len(keys) != len(set(keys)) or not rows:
        raise ValueError("MuData schema catalog is empty or has duplicate fields")
    print(f"Prepared {len(PAGES) + 2} documentation pages and {len(rows)} MuData fields")


if __name__ == "__main__":
    main()
