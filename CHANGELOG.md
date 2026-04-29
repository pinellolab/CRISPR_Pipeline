# nf-core/crispr: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - 2026-April-29

### `Added`

- Added FASTQ subsampling parameters: `subsample_seed`, `subsample_probability`, and `subsample_record_count`.
- Added `conf/modules/subsample.config` to configure `FQ_SUBSAMPLE` arguments, output prefix, and publish directory.
- Added local `FASTQ_SUBSAMPLE` subworkflow wrapping the nf-core `fq/subsample` module and emitting subsampled reads plus module versions.

## v1.0.0 - 2025-AUG-19

- Initial release of IGVF Perturb-seq Pipeline
- Core pipeline functionality
- Documentation and examples

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
