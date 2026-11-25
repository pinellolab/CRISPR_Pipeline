# Documentation sources

Read the [CRISPR Pipeline website](https://pinellolab.github.io/CRISPR_Pipeline/)
for usage, current dataset examples, analysis, outputs, and parameters.
GitHub Actions builds and deploys the site from `main`.

To preview the website locally, run:

```bash
python3 scripts/check_docs_sync.py
python3 scripts/build_docs_site.py
mkdocs build --strict
```

Canonical prose is in this directory. The site downloads are built from the
versioned configs and `examples/datasets/` in the repository. The parameter
index comes from `nextflow_schema.json`, and the MuData field dictionary comes
from `mudata_schema_catalog.tsv`. See the repository documentation-sync skill
for checks and the update checklist.
