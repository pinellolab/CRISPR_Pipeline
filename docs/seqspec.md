# SeqSpec for RNA, guides, and HTO

SeqSpec describes where a library stores the cell barcode, UMI, transcript/guide/HTO feature, and any fixed sequence. Supply the correct YAML path in each samplesheet row; the pipeline validates guide and optional HTO SeqSpecs before mapping. **Do not copy the example YAML unchanged for a different chemistry.**

| Modality | Example | What to inspect |
|---|---|---|
| RNA | [RNA SeqSpec](downloadable-examples/rna_seqspec.yml) | `sequence_spec` R1/R2 identifiers, barcode and UMI lengths/positions, cDNA orientation, onlist. |
| Guide | [Guide SeqSpec](downloadable-examples/guide_seqspec.yml) | Guide read ID and strand, barcode/UMI, captured spacer region, fixed adapters. |
| HTO | [Hash SeqSpec](downloadable-examples/hash_seqspec.yml) | HTO read and feature region, barcode/UMI, HTO reference. Required only when hashing is enabled. |

The repository examples show 16-base barcodes and 12-base UMIs for one historical assay. Those numbers are **not** universal defaults. Compare the YAML to a sample of the actual reads and to the sequencing provider's library design. A spacer tag configured through `spacer_tag` asks guide mapping to search the full guide read around that fixed sequence; it does not repair an incorrect read orientation, barcode definition, or guide reference. `reverse_complement_guides` is a separate capture-orientation choice and should be supported by read evidence.

For an assay with a few measurement sets, test one complete set using `DEMO_MODE=true`, then inspect the SeqSpec check, guide-mapping recovery, and barcode-overlap sections of the dashboard before running all sets. Guide metadata sequences must describe the **captured** guide sequence after any justified orientation/capture transformation, with the original biological spacer retained in a separate provenance column when applicable.
