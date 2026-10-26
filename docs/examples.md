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

## Completed benchmark output example

The [Gary Hon full-depth benchmark example on Synapse](https://www.synapse.org/Synapse:syn78349990)
shows an actual completed pipeline run (`garyhondataset_dashboard_calibration_20261022_01`,
source commit `4aa7cc6`). It includes the full inference MuData, local/global
per-guide and per-element Parquet tables, a bundled QC archive, the final
dashboard, inputs, and provenance. Raw FASTQs and Nextflow work files are not
duplicated there. Use the included README and manifest to interpret the
package; its samplesheet records the original server paths and must be
relocalized before reuse.
