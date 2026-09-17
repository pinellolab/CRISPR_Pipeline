# Per-measurement-set RNA QC

RNA cell calling and cell-level QC run independently for every samplesheet
`measurement_sets` value before the retained matrices are concatenated. This
prevents a deep lane from determining the knee or robust thresholds for a
shallower lane. The independent Nextflow tasks can run in parallel and use the
pipeline base container; the workflow does not create a host Python environment.

For each measurement set the pipeline:

1. reads its mapped unfiltered RNA AnnData;
2. qualifies cell barcodes with the same measurement-set key used by the guide
   and hashing modalities;
3. computes and plots its own barcode-rank curve and `knee`/`knee2` points;
4. applies the selected knee, fixed RNA UMI/minimum-gene/mitochondrial filters,
   and any enabled MAD filters;
5. writes one filtered AnnData, two QC plots, and one audit-table row; and
6. concatenates only the retained measurement-set matrices.

The standard 10-cell gene-support floor (or one-cell floor in TAP-seq mode) is
applied after concatenation. Gene support remains global intentionally: a real
gene is not removed merely because it is sparse in one measurement set. The
later fractional `QC_min_cells_per_gene` filter is unchanged.

## MAD options

All MAD options are disabled by default (`0`) and can be enabled independently:

| Parameter | Metric | Tail |
|---|---|---|
| `QC_MAD_total_counts` | `log1p(total_counts)` | lower and upper |
| `QC_MAD_n_genes` | `log1p(n_genes_by_counts)` | lower and upper |
| `QC_MAD_pct_mito` | `pct_counts_mt` | upper only |

For a value `k`, a two-sided metric retains values in
`median - k*MAD` through `median + k*MAD`; mitochondrial filtering has no MAD
lower bound. A metric with zero MAD is not filtered, avoiding removal of an
entire tied population. MAD filters are combined with the fixed thresholds, so
enabling them never relaxes `QC_min_counts_per_cell`, `QC_min_genes_per_cell`,
or `QC_pct_mito`.

Example:

```groovy
params {
    QC_MAD_total_counts = 3
    QC_MAD_n_genes = 3
    QC_MAD_pct_mito = 3
}
```

## Outputs

The combined artifacts are published immediately under
`measurement_set_qc/` so the live W&B execution dashboard can display them
before the final dashboard process runs. The final dashboard `figures/`
directory also contains:

- `knee_plot_scRNA_<measurement_set>.png`;
- `qc_distributions_scRNA_<measurement_set>.png`;
- `rna_qc_filter_flow_<measurement_set>.png`, showing cells → resolved filter
  parameter → retained cells in exact execution order;
- `rna_qc_filter_steps_<measurement_set>.png`, showing the relevant RNA metric
  immediately before and after every fixed or MAD filter;
- `measurement_set_qc_filter_flow.tsv`, the machine-readable sequential counts,
  thresholds, enabled/skipped state, removal percentage, and retained percentage; and
- `measurement_set_qc_metrics.tsv`, with input, post-knee, post-fixed-threshold,
  retained-cell, threshold, median, MAD, and bound values for every measurement
  set.

The dashboard RNA-QC block renders these plots and the combined table.
The displayed order is barcode calling (`none`, `knee`, or `knee2`), minimum
RNA UMIs, minimum detected genes (only for `none`), fixed mitochondrial cutoff,
RNA-UMI MAD, detected-gene MAD, and mitochondrial MAD. Disabled or inapplicable
steps remain visible with zero removal. The global gene-support filter remains
separate because it runs after measurement-set concatenation.
