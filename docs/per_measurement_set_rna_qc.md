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
4. applies the selected knee, the fixed RNA UMI floor, and two-sided RNA UMI
   and detected-gene MAD filters;
5. optionally runs Scrublet and removes predicted doublets independently within
   that measurement set;
6. concatenates retained measurement sets, then applies the fixed mitochondrial
   cell cutoff and fractional gene-support filter;
7. writes filtered AnnData objects, QC plots, and machine-readable audit tables.

Gene support remains global intentionally: a real gene is not removed merely
because it is sparse in one measurement set. The only gene-support rule is the
post-concatenation fractional `QC_min_cells_per_gene` threshold.

## MAD options

The RNA-complexity MAD filters default to `5` and can be changed independently:

| Parameter | Metric | Tail |
|---|---|---|
| `QC_MAD_total_counts` | `log1p(total_counts)` | lower and upper |
| `QC_MAD_n_genes` | `log1p(n_genes_by_counts)` | lower and upper |

For a value `k`, a two-sided metric retains values in
`median - k*MAD` through `median + k*MAD`. A metric with zero observed MAD is
not filtered, avoiding removal of an entire tied population. Mitochondrial MAD
filtering and the fixed per-cell detected-gene threshold are no longer active.

Example:

```groovy
params {
    QC_MAD_total_counts = 5
    QC_MAD_n_genes = 5
    ENABLE_SCRUBLET = true
    SCRUBLET_assay_type = 'droplet'
    SCRUBLET_expected_doublet_rate = null
    SCRUBLET_n_prin_comps = 30
    SCRUBLET_adaptive_pca_fallback = true
}
```

## Outputs

The combined artifacts are published immediately under
`measurement_set_qc/` so the live W&B execution dashboard can display them
before the final dashboard process runs. The final dashboard `figures/`
directory also contains:

- `knee_plot_scRNA_<measurement_set>.png`;
- `qc_distributions_scRNA_<measurement_set>.png`;
- `rna_qc_filter_flow_<measurement_set>.png`, with one separate
  cells-before → filter → cells-after row per step in exact execution order;
- `rna_qc_filter_steps_<measurement_set>.png`, whose first row shows the full
  barcode-rank curve, only the knee actually applied, the retained curve, and
  explicit before/after cell counts. Later rows show before/after histograms
  and boxplots for each UMI, MAD, or Scrublet filter. Two-sided filters label
  their lower and upper MAD boundaries;
- `scrublet_scores_scRNA_<measurement_set>.png` when Scrublet is enabled;
- `post_concat_mito_before_after.png`, `post_concat_gene_support.png`, and
  `post_concat_qc_filter_flow.png` for the filters applied after concatenation;
- `measurement_set_qc_filter_flow.tsv`, the machine-readable sequential counts,
  thresholds, enabled/skipped state, removal percentage, and retained percentage; and
- `measurement_set_qc_metrics.tsv`, with input, post-knee, post-fixed-threshold,
  retained-cell, threshold, median, MAD, and bound values for every measurement
  set.

The dashboard RNA-QC block renders these plots and the combined table.
The displayed order is barcode calling (`none`, `knee`, or `knee2`), minimum
RNA UMIs, RNA-UMI MAD, detected-gene MAD, and Scrublet. After measurement-set
concatenation, `QC_pct_mito` is applied to cells and then
`QC_min_cells_per_gene` is applied to genes. Disabled or inapplicable steps
remain visible with zero removal.

`SCRUBLET_assay_type = 'droplet'` resolves an automatic expected-doublet rate
of `0.08`; `cc-perturb-seq` resolves to `0.025`. A numeric
`SCRUBLET_expected_doublet_rate` overrides the profile.

When `SCRUBLET_adaptive_pca_fallback=true`, only Scrublet's explicit
`n_components` dimensionality error is retried. The retry uses one fewer
component than Scrublet's reported usable dimension. Requested and actual PCA
dimensions and whether fallback occurred are saved in
`measurement_set_qc_metrics.tsv`; unrelated Scrublet errors remain fatal.
