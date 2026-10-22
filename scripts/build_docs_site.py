#!/usr/bin/env python3
"""Assemble the GitHub Pages source from canonical checked-in documentation.

No pipeline environment or network connection is needed for this step.
"""

import csv
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
DEST = ROOT / "build" / "docs-site"

PAGES = (
    "usage.md",
    "samplesheet.md",
    "seqspec.md",
    "guide-metadata.md",
    "samplesheet_covariates.md",
    "run-recipes.md",
    "preinference_checkpoint_user_guide.md",
    "filtering_workflow.md",
    "postconcat_embedding_qc.md",
    "inference.md",
    "live_wandb_execution_dashboard.md",
    "output.md",
    "mudata_schema.md",
    "qc_outputs.md",
    "target_aware_gene_filter.md",
    "guide_mapping_overlap_qc.md",
    "clone_and_saturation_qc.md",
    "examples.md",
)


def flatten_properties(node, prefix=""):
    for key, value in node.get("properties", {}).items():
        name = f"{prefix}.{key}" if prefix else key
        if "properties" in value:
            yield from flatten_properties(value, name)
        else:
            yield name, value


def render_parameters():
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    lines = [
        "# Parameter reference",
        "",
        "Generated from `nextflow_schema.json` by `scripts/build_docs_site.py`.",
        "Values can also be supplied in a Nextflow config or a JSON/YAML params file.",
        "Dataset-specific values must be reviewed before launch.",
        "",
        "| Parameter | Type | Schema default | Description |",
        "|---|---|---|---|",
    ]
    parameters = {}
    for section in schema.get("$defs", {}).values():
        for name, spec in flatten_properties(section):
            parameters.setdefault(name, spec)
    if not parameters:
        raise ValueError("No parameters found in nextflow_schema.json")
    for name, spec in parameters.items():
        value = spec.get("default", "—")
        default = json.dumps(value, ensure_ascii=False) if value != "—" else "—"
        description = spec.get("description", "").replace("|", "\\|").replace("\n", " ")
        help_text = spec.get("help_text", "").replace("|", "\\|").replace("\n", " ")
        desc = description + (f" {help_text}" if help_text else "")
        lines.append(f"| `{name}` | `{spec.get('type', '—')}` | `{default}` | {desc} |")
    lines.insert(5, f"This index contains **{len(parameters)}** schema parameters. ")
    return "\n".join(lines) + "\n"


def main():
    DEST.mkdir(parents=True, exist_ok=True)
    for name in PAGES:
        shutil.copyfile(DOCS / name, DEST / name)
    shutil.copyfile(DOCS / "website" / "index.md", DEST / "index.md")
    shutil.copyfile(DOCS / "mudata_schema_catalog.tsv", DEST / "mudata_schema_catalog.tsv")
    (DEST / "parameters.md").write_text(render_parameters())

    examples = DEST / "downloadable-examples"
    examples.mkdir(exist_ok=True)
    for source, target in (
        ("docs/examples/minimal_samplesheet.csv", "samplesheet.csv"),
        ("docs/examples/hashing_samplesheet.csv", "hashing_samplesheet.csv"),
        ("example_data/guide_metadata.tsv", "guide_metadata.tsv"),
        ("example_data/hash_metadata.tsv", "hash_metadata.tsv"),
        ("example_data/yaml_files/rna_seqspec.yml", "rna_seqspec.yml"),
        ("example_data/yaml_files/guide_seqspec.yml", "guide_seqspec.yml"),
        ("example_data/yaml_files/hash_seqspec.yml", "hash_seqspec.yml"),
    ):
        shutil.copyfile(ROOT / source, examples / target)
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
