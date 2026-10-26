---
name: sync-pipeline-docs
description: Keep CRISPR Pipeline source, parameter and MuData contracts, plugin installation guidance, and the GitHub Pages documentation synchronized when changing pipeline behavior, public inputs or outputs, dashboards, or documentation.
---

# Sync pipeline and documentation

Use this skill for any change to public parameters/defaults, samplesheet or metadata
contracts, QC/inference behavior, output fields, dashboard/plugin behavior, run
instructions, or GitHub Pages content. It is a review workflow, not a guarantee
that prose is automatically correct.

1. Read the affected implementation and tests. Identify the actual source of
   truth: `nextflow.config` and `nextflow_schema.json` for parameters;
   workflow/modules/scripts for behavior; `docs/mudata_schema_catalog.tsv`
   for published MuData fields; `conf/dashboard_plugin.config`,
   `run_support/install_dashboard_plugin.sh`, and the bundled ZIP for the
   plugin contract.
2. Update the relevant canonical `docs/*.md` page and examples in the same
   change. For new pages, add them to both `mkdocs.yml` and
   `scripts/build_docs_site.py`. Explain any default changes, migration or
   resume implications, and assay-specific exceptions. Keep the website and
   repository README linked to the canonical page; avoid copy-pasted competing
   instructions.
3. For plugin changes, verify version and ZIP digest against the installer and
   config. Keep `docs/plugin-installation.md` accurate. Do not say a fresh
   clone installs the plugin or that registry publication exists without
   testing those claims. Do not embed secrets or user-specific paths.
4. Run `python scripts/check_docs_sync.py`,
   `python scripts/build_docs_site.py`, and `mkdocs build --strict` (in an
   environment with MkDocs installed). Run focused tests, including
   `tests/test_docs_site.py` and `tests/test_docs_sync.py`. Review the rendered
   page and links. If changing the MuData catalog, regenerate its rendered
   reference with `bin/render_mudata_schema_docs.py` and inspect the diff.
5. Before a commit or PR, review code, schema, docs, site, tests, and plugin
   installation instructions together. CI catches mechanical drift; explicitly
   report any semantic behavior that still needs a manual or end-to-end check.
   Do not push, publish, or alter a live run unless the user requested it.

On the Pinello server, preserve clean pinned pipeline checkouts and follow
`AGENTS.md` local-run storage and provenance rules. This skill never authorizes
deleting work directories, caches, or running jobs.
