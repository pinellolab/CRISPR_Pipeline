# Gasperini at-scale K562 CRISPRi

For CRISPR Pipeline **v1.1.0**. Based on the full-input November 2026
`rerun_paper_filters_20261123` analysis. The samplesheet covers 32 measurement
sets and 388 RNA/guide rows. Guide coordinates are in hg38.

- Raw reads and study accession: [GSE120861](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861).
- `sequencing_run` identifies the original SRA record. The input FASTQ names
  retain the relative BAM-to-FASTQ output path and lane. Recover the
  barcode/UMI and cDNA or guide read from the original BAM with
  [10x bamtofastq](https://github.com/10XGenomics/bamtofastq), and match each
  converted lane to the paths and read roles in the samplesheet.
- Supply `737K-august-2016.txt` from the Cell Ranger 3′ v2 whitelist.
  Expected uncompressed MD5 is `f62a276e262fdd85262a889d0f48556b`.
- The config carries over the RNA and guide-assignment settings from the
  source run. The paper-specific hypothesis set and external evaluation
  are separate. Full mode uses the release's standard inference settings.

Use the localization and launch commands in the
[examples guide](https://pinellolab.github.io/CRISPR_Pipeline/examples/).
