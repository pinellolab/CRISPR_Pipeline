#!/usr/bin/env python3
"""Fail when the published site or bundled plugin instructions drift from source."""

import hashlib
import json
import re
from pathlib import Path

from build_docs_site import DOCS, PAGES, ROOT, render_configuration
from render_slurm_config import OUTPUT, render as render_slurm_config


def require(condition, message):
    if not condition:
        raise SystemExit(message)


def capture(pattern, source, label):
    match = re.search(pattern, source, re.MULTILINE)
    require(match is not None, f"Cannot locate {label}")
    return match.group(1)


def main():
    nav = (ROOT / "mkdocs.yml").read_text()
    listed = set(re.findall(r"(?<![\w/-])([\w-]+\.md)\s*$", nav, re.MULTILINE))
    expected = set(PAGES) | {"index.md"}
    require(listed == expected, f"Site navigation differs from copied pages: missing={expected-listed}, extra={listed-expected}")
    require(all((DOCS / name).is_file() for name in PAGES), "A site source page is missing")
    configuration = render_configuration()
    config_text = (ROOT / "nextflow.config").read_text()
    gate = "GUIDE_MAPPING_QC_min_overlap_to_rna_fraction"
    configured = capture(rf"^\s*{gate}\s*=\s*([0-9.]+)$", config_text, "guide overlap default")
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    schema_defaults = [spec["default"] for section in schema["$defs"].values()
                       for name, spec in section.get("properties", {}).items() if name == gate]
    require(len(schema_defaults) == 1 and float(configured) == schema_defaults[0], "Guide overlap default differs from schema")
    require(f"defaults to {configured}" in (DOCS / "guide_mapping_overlap_qc.md").read_text(), "Guide overlap guide has stale default")
    for name, page in (("ANALYSIS_mode", "usage.md"),
                       ("QC_barcode_filter", "configuration.md"),
                       ("INFERENCE_PERTURBO_GLOBAL_RESULTS_FORMAT", "output.md")):
        configured_value = capture(rf"^\s*{name}\s*=\s*'([^']+)'$", config_text, f"{name} default")
        values = [spec["default"] for section in schema["$defs"].values()
                  for key, spec in section.get("properties", {}).items() if key == name]
        require(len(values) == 1 and configured_value == values[0], f"{name} default differs from schema")
        content = configuration if page == "configuration.md" else (DOCS / page).read_text()
        require(configured_value in content, f"{page} does not mention {name} default")

    installer = (ROOT / "run_support/install_dashboard_plugin.sh").read_text()
    config = (ROOT / "conf/dashboard_plugin.config").read_text()
    build = (ROOT / "plugins/nf-crispr-dashboard/build.gradle").read_text()
    guide = (DOCS / "plugin-installation.md").read_text()
    version = capture(r"^VERSION=([0-9.]+)$", installer, "installer version")
    config_version = capture(r"id 'nf-crispr-dashboard@([0-9.]+)'", config, "config version")
    build_version = capture(r"^version = '([0-9.]+)'$", build, "build version")
    require(version == config_version == build_version, "Plugin versions differ between installer, config, and build")
    archive = ROOT / f"plugins/nf-crispr-dashboard/dist/nf-crispr-dashboard-{version}.zip"
    require(archive.is_file(), f"Bundled plugin archive missing: {archive}")
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    pinned = capture(r"^EXPECTED_SHA256=([0-9a-f]{64})$", installer, "installer checksum")
    require(digest == pinned, "Bundled plugin archive differs from pinned installer checksum")
    require(all(value in guide for value in (version, digest, "NXF_PLUGINS_DIR", "run_support/install_dashboard_plugin.sh", "conf/dashboard_plugin.config")), "Plugin installation guide is out of sync")
    require("plugin-installation.md" in (DOCS / "website/index.md").read_text(), "Website landing page does not link plugin installation")
    require(OUTPUT.is_file() and OUTPUT.read_text() == render_slurm_config(),
            "Slurm defaults are stale. Run python scripts/render_slurm_config.py")
    manifest = (ROOT / "conf/provenance.config").read_text()
    release = capture(r"^\s*version = '([^']+)'", manifest, "pipeline manifest version")
    for relative in ("nextflow.config", "nextflow_cc.config", "nextflow_tapseq.config",
                     "conf/simplified_pipeline/nextflow_yang.dev.config", "docs/examples.md"):
        require(f"v{release}" in (ROOT / relative).read_text(), f"{relative} has a stale release version")
    print(f"Documentation and bundled plugin checks passed ({len(PAGES)} pages, plugin {version})")


if __name__ == "__main__":
    main()
