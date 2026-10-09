# RNA–guide barcode overlap QC

The default recovery gate uses `overlap_cells / rna_cells`: the fraction of
high-quality RNA cells represented in the guide matrix, calculated separately
for each measurement set. `GUIDE_MAPPING_QC_min_overlap_to_rna_fraction`
defaults to 0.5. Minimum overlapping cells (20) and expected guide-feature
recovery (0.5) remain independently enforced.

RNA has already undergone cell calling and RNA QC at this stage. The raw guide
matrix may still contain many uncalled droplets. Dividing the shared barcodes
by all raw guide barcodes penalizes valid RNA filtering and is therefore not
the default gate. `overlap_to_guide_fraction` and raw guide/RNA counts are still
published for diagnosis. The optional legacy threshold
`GUIDE_MAPPING_QC_min_overlap_to_guide_fraction` defaults to 0; explicitly
setting a nonzero value retains its previous raw-guide-denominator behavior.

This is a barcode-presence QC, not a guide-assignment rate. Presence in the
mapped matrix does not establish a positive guide assignment. Guide assignment
and its cell filters run afterwards; low-overlap or insufficient guide-feature
recovery is still an error by default.
