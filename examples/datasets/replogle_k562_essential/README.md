# Replogle K562 day-6 essential-gene Perturb-seq

For CRISPR Pipeline **v1.1.0**. The samplesheet comes from the October 2026
pre-inference run that supplied the subsequent full inference.
It covers 48 GEM groups with 384 RNA and 192 guide rows. The included
hg38 guide metadata preserves the dual-guide construct identities.

- Raw reads: [SRA BioProject PRJNA831566](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA831566).
  Use the `sequencing_run` accessions and the paths in the samplesheet to
  stage each mate under the correct library folder. Example names retain
  the `SRR…_1.fastq.gz` and `SRR…_2.fastq.gz` suffixes.
- Study portal: [K562 essential](https://gwps.wi.mit.edu/K562_essential).
  The input list covers day 6, with separately identified RNA and sgRNA runs.
- Supply the uncompressed 10x Chromium 3′ v3 `3M-february-2018.txt`
  barcode allowlist. Expected MD5 is `3d36d0a4021fd292b265e2b5e72aaaf3`.

Use the localization and launch commands in the
[examples guide](https://pinellolab.github.io/CRISPR_Pipeline/examples/).
