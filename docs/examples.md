# Downloadable input examples

These are format references checked into the repository, **not** a ready-to-run dataset. The samplesheet has illustrative absolute paths and the SeqSpecs/metadata describe a specific assay. Replace paths and verify library structure, guide orientation, barcode chemistry, reference build, and captured sequences for your experiment.

| File | Purpose |
|---|---|
| [Samplesheet](downloadable-examples/samplesheet.csv) | Current required columns for one RNA/guide measurement set. |
| [HTO samplesheet](downloadable-examples/hashing_samplesheet.csv) | Adds an HTO row and metadata map for one set. |
| [RNA SeqSpec](downloadable-examples/rna_seqspec.yml) | Cell barcode, UMI and transcript structure. |
| [Guide SeqSpec](downloadable-examples/guide_seqspec.yml) | Guide read structure. |
| [HTO SeqSpec](downloadable-examples/hash_seqspec.yml) | Optional hashing read structure. |
| [Guide metadata](downloadable-examples/guide_metadata.tsv) | Guide IDs, spacers and intended targets. |
| [HTO metadata](downloadable-examples/hash_metadata.tsv) | HTO feature reference. |

See [run recipes](run-recipes.md) for command examples and [all parameters](parameters.md) for the current configurable surface.
