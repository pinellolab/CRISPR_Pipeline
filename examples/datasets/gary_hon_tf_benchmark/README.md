# Gary Hon WTC11 TF perturbation benchmark

For CRISPR Pipeline **v1.1.0**. The input and settings come from the completed
October 2026 full-depth TF benchmark `garyhondataset_dashboard_calibration_20261022_01`.
The samplesheet covers three measurement sets with RNA, guide, and hashtag
libraries. Separate HTO metadata maps the captured hashtags to samples.

- [IGVF dataset measurement set IGVFDS9613DDRB](https://data.igvf.org/measurement-sets/IGVFDS9613DDRB/)
  and two additional sets listed in `sample_metadata.tsv`.
- `fastq_accessions.csv` gives IGVF identifiers, download URLs, sizes and
  MD5 checksums. Keep the accession-based FASTQ filenames.
- Supply barcode allowlist `IGVFFI9487JPEN.tsv` from the IGVF record.
  Expected MD5 is `f62a276e262fdd85262a889d0f48556b`.
- [Completed benchmark package](https://www.synapse.org/Synapse:syn78349990).

Use the localization and launch commands in the
[examples guide](https://pinellolab.github.io/CRISPR_Pipeline/examples/).
