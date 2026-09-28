# Target-aware gene retention

The global RNA prevalence filter remains the default. Set
`QC_target_gene_rescue = true` to retain an additional low-prevalence gene only
when all three conditions hold:

```text
canonical intended target
AND unique assigned target-positive cells >= QC_target_gene_min_assigned_cells
AND detected RNA fraction >= QC_target_gene_min_detected_fraction
```

The full rule is the global `QC_min_cells_per_gene` condition OR this rescue.
Detection is calculated from raw RNA counts after final cell-level QC. Guide
support is calculated from binary `guide.layers['guide_assignment']`, with
cells deduplicated across guides targeting the same gene. Non-targeting,
safe-targeting, and negative-control guides cannot rescue genes.

Defaults preserve historical behavior:

- `QC_target_gene_rescue = false`
- `QC_target_gene_min_assigned_cells = 20`
- `QC_target_gene_min_detected_fraction = 0.001` (0.1%)

The post-concatenation QC output includes `gene_filter_decisions.tsv.gz` with
the prevalence, assignment support, decision, and reason for every input gene.
Retained genes carry the same decision fields in `gene.var`.
