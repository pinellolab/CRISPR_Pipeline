# CRISPR Pipeline documentation

The [GitHub Pages website](https://pinellolab.github.io/CRISPR_Pipeline/) is built from these checked-in pages and the pipeline schema. Its entry points are [start here](usage.md), [samplesheet](samplesheet.md), [SeqSpec](seqspec.md), [guide metadata](guide-metadata.md), [run recipes](run-recipes.md), and [pause/edit/continue](preinference_checkpoint_user_guide.md). The complete [MuData field dictionary](mudata_schema.md) is generated from `mudata_schema_catalog.tsv`; the website's parameter index is generated from `nextflow_schema.json`.

To build the website locally: `python3 scripts/build_docs_site.py && mkdocs build --strict`. GitHub Actions runs the same checks on `dev` and deploys to Pages. Do not edit generated `build/docs-site` or `build/site` files.

The nf-core/crispr documentation is split into the following pages:

- [Usage](usage.md)
  - An overview of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.
- [Filtering and MuData assembly](filtering_workflow.md)
  - The filtering, barcode intersection, guide assignment, and concatenation order for multiple measurement sets.
- [MuData field reference](mudata_schema.md)
  - Field-level data dictionary for the final MuData object, including modalities, cell and feature annotations, assignment layers, optional hashing fields, and inference result tables.
- [Interactive filtering workflow](filtering_workflow.html)
  - Standalone HTML with rendered Mermaid diagrams, implementation snippets, and the active CC configuration.
- [Parameter reference](index.html)
  - Standalone nf-core-style parameter reference generated from `nextflow_schema.json`.

Regenerate the parameter reference with:

```bash
python3 bin/render_schema_index.py
```

The generator always writes `docs/index.html` so the published page keeps the same filename.

You can find a lot more documentation about installing, configuring and running nf-core pipelines on the website: [https://nf-co.re](https://nf-co.re)
