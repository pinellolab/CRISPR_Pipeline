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

## Fixed guide intervals in FLASH base-editing mapping

`BASEEDITING_FLASH_search_mode = 'full_read'` retains the legacy sliding search.
For assays whose SeqSpec defines the captured spacer exactly, select `seqspec`:
only that interval is compared with guide references, and reads tied between
multiple best-matching guides are left unassigned. The interval length must equal
the matched guide length, and captured reference sequences must be unique.
The mode and rejected ambiguous-read count are recorded in the mapped AnnData
parameters. Alignment QC positions in this mode are relative to the extracted
interval. Barcode and UMI extraction remain defined by the same SeqSpec.

A full-read sliding search can find overlapping, shifted spacers from different
guide designs at different positions in one read. Even zero-mismatch matches
can therefore change guide identity. Use the assay's validated feature interval
when comparing such libraries with the original publication; lowering mismatch
tolerance alone does not resolve cross-position ambiguity. Changing this mode
requires guide remapping and downstream assignment/inference, rather than only
regenerating evaluation plots. This option is unrelated to RNA cell calling.
