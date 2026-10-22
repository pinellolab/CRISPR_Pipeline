#!/usr/bin/env python3
"""Fail when the published site or bundled plugin instructions drift from source."""

import hashlib
import json
import re
from pathlib import Path

from build_docs_site import DOCS, PAGES, ROOT, render_parameters


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
    expected = set(PAGES) | {"index.md", "parameters.md"}
    require(listed == expected, f"Site navigation differs from copied pages: missing={expected-listed}, extra={listed-expected}")
    require(all((DOCS / name).is_file() for name in PAGES), "A site source page is missing")
    require("| `ANALYSIS_mode` |" in render_parameters(), "Generated parameter reference lacks ANALYSIS_mode")
    gate = "GUIDE_MAPPING_QC_min_overlap_to_rna_fraction"
    configured = capture(rf"^\s*{gate}\s*=\s*([0-9.]+)$", (ROOT / "nextflow.config").read_text(), "guide overlap default")
    schema = json.loads((ROOT / "nextflow_schema.json").read_text())
    schema_defaults = [spec["default"] for section in schema["$defs"].values()
                       for name, spec in section.get("properties", {}).items() if name == gate]
    require(len(schema_defaults) == 1 and float(configured) == schema_defaults[0], "Guide overlap default differs from schema")
    require(f"defaults to {configured}" in (DOCS / "guide_mapping_overlap_qc.md").read_text(), "Guide overlap guide has stale default")

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
    print(f"Documentation and bundled plugin checks passed ({len(PAGES)} pages, plugin {version})")


if __name__ == "__main__":
    main()
