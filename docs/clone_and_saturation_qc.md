# Clone filtering and sequencing-saturation QC

Both analyses run in the pipeline base container. They are disabled by default
and can be enabled independently.

## Guide-barcode clone detection

Enable clone detection after guide assignment and before inference:

```groovy
params {
    ENABLE_CLONE_REMOVAL = true
    CLONE_REMOVAL_action = 'drop_clonal'
    CLONE_REMOVAL_alpha = 0.05
    CLONE_REMOVAL_min_clone_size = 2
}
```

The implementation follows the agglomerative hypergeometric method described
by Wang et al. and its MIT-licensed reference implementation:

- Paper: https://doi.org/10.1186/s12864-022-08359-1
- Reference scripts: https://github.com/yihan1119/Group_clone/tree/main/Scripts
- Copied license: `third_party/Group_clone_LICENSE`

For every cell, the set of positive values in
`guide.layers['guide_assignment']` is treated as its guide barcode. The cell is
compared with the representative cell of each previously identified clone.
The overlap is tested with `scipy.stats.hypergeom.sf(overlap - 1, M, n, N)`
and the alpha is Bonferroni-corrected over all possible cell pairs. A cell
matching no clone starts a clone, a cell matching one clone joins it, and a
cell matching multiple clones is labeled ambiguous/doublet.

Like the reference implementation, grouping is sequential and therefore uses
the first cell of each discovered group as its representative. The pipeline
makes the traversal deterministic by preserving the MuData cell order and
writes every assignment to an audit table.

Actions are:

- `mark_only`: preserve every cell and annotate the retained MuData.
- `drop_clonal`: remove every member of clone groups at least
  `CLONE_REMOVAL_min_clone_size`, plus ambiguous cells.
- `keep_representative`: retain the cell with the largest RNA UMI total in
  each detected clone and remove its other members plus ambiguous cells.

This method is not a generic low-MOI duplicate detector. The source paper used
approximately 20,000 guides and a median of 32 guides per cell and recommends
at least 1,000 guides and at least 10 guides per cell. The pipeline therefore
keeps it opt-in and reports `low_power_warning` in the dashboard when that rule
of thumb is not met.

Outputs are published under `clone_removal/` and copied into the dashboard:

- `clone_cell_assignments.tsv`: every input cell, clone, status, and retention.
- `clone_groups.tsv`: clone sizes and retained counts.
- `clone_metrics.tsv` and `.json`: run-level filtering and applicability QC.
- `clone_filter_summary.png`: clone-size and filter-outcome visualization.

## 10x-style sequencing saturation

Enable RNA saturation curves:

```groovy
params {
    ENABLE_SEQUENCING_SATURATION = true
    SATURATION_downsample_fractions = '0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0'
}
```

The implementation reads each mapping batch's corrected and sorted
`output.unfiltered.bus`. It restricts the calculation to cell barcodes retained
by RNA preprocessing and to equivalence classes mapping unambiguously to one
gene. For a molecule supported by `r` usable reads, its probability of being
observed at read fraction `p` is calculated analytically as
`1 - (1 - p)^r`. No downsampled FASTQs are written.

At each depth, sequencing saturation is:

```text
1 - expected unique (cell barcode, UMI, gene) molecules / expected usable reads
```

This matches the duplicate-read definition documented by 10x Genomics. The
curves also report median expected UMIs and genes per retained cell. These are
RNA mapping curves; the pipeline does not substitute RNA saturation for guide
capture saturation.

Outputs are published under `sequencing_saturation/` and copied into the
dashboard:

- `sequencing_saturation_curve.tsv`: aggregate and per-batch rarefaction data.
- `sequencing_saturation_metrics.tsv`: full-depth endpoint metrics.
- `sequencing_saturation_curve.png`: saturation, UMI, and gene curves.
- `sequencing_saturation_method.json`: exact population and mapping rules.

10x definitions:

- https://www.10xgenomics.com/support/software/cell-ranger/9.0/analysis/outputs/cr-3p-outputs-metrics-count
- https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/web-summary
