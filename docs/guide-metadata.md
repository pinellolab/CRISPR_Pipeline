# Guide and HTO metadata

The guide metadata connects each captured guide sequence to a stable identifier, targeting status, and (for cis/local analysis) its intended genomic element. The [repository example](downloadable-examples/guide_metadata.tsv) is a TSV with these columns:

```text
guide_id  spacer  targeting  type  guide_chr  guide_start  guide_end  strand  pam  intended_target_name  intended_target_chr  intended_target_start  intended_target_end
```

`guide_id` must be unique. `spacer` must match the sequence the mapping step can see, not merely a design-sheet alias. `targeting` should be a real boolean and `type` should distinguish targeting from non-targeting/control guides. Preserve controls: they are needed for some guide-calling checks, low-MOI contrasts, and non-targeting calibration. For a targeting guide, provide consistent `intended_target_name` and element coordinates so guides for the same element aggregate correctly. A guide's genomic coordinates and its intended element coordinates are **different concepts**; do not fill one with the other. For controls, avoid invented genomic positions.

Coordinate values must be in one declared reference build and convention compatible with the GTF used for inference. Record the source build and any liftover or metadata repair in additional columns; the pipeline preserves source-specific guide columns in `guide.var`. Gene symbols and Ensembl IDs need separate columns when both are available. Never change source metadata merely to satisfy a downstream join without documenting the transformation.

Validation checklist:

1. No duplicate guide IDs or ambiguous captured sequences unless the assay intentionally allows them.
2. Targeting and control classifications agree with design records; non-targeting controls remain present.
3. Captured spacer length, orientation and any tag are confirmed from actual guide reads.
4. Intended-target names are stable within an element; coordinates, chromosome prefix and GTF build agree.
5. Negative controls do not acquire accidental intended-target genes or coordinates.
6. For dual-guide assays, retain both guide identities and enable `DUAL_GUIDE`; expect the post-assignment collapse to select the intended two-guide cells.

For hashing assays, the `barcode_hashtag_map` file supplies HTO feature IDs/sequences (see the [example](downloadable-examples/hash_metadata.tsv)). The HTO read SeqSpec, reference, and mapping must match each other. HTO singlet and minimum-positive-cell filters are applied to the already RNA/guide-filtered population; a weak or mis-oriented HTO library can remove many usable-looking RNA cells. Compare per-set HTO counts and GMM-Demux calls in the dashboard before relaxing thresholds.
