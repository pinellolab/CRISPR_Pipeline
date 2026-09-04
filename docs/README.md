# nf-core/crispr: Documentation

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
