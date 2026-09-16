# nf-core/crispr: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Changed`

- **One control group for both inference methods.** The new `INFERENCE_control_group`
  parameter (default `auto`) names the cells a perturbation is compared against,
  in SCEPTRE's vocabulary, and the pipeline translates it for PerTurbo:
  `nt_cells` -> `--crt-pool control-anchored`, `complement` -> `--crt-pool all-cells`,
  `auto` -> whichever `Multiplicity_of_infection` implies. The resolution is logged
  once at the start of the inference subworkflow and recorded beside the results
  (`control_group_resolution.json`, a `pipeline_control_group` block in PerTurbo's
  `crt_metadata.json`, and `sceptre_control_group.json`).

### `Fixed`

- **SCEPTRE now honours its control group, so SCEPTRE results on low-MOI screens
  change.** `INFERENCE_SCEPTRE_control_group` was threaded into the SCEPTRE driver
  and then discarded: the driver always used the complement contrast and emitted
  only a warning. Under the `auto` default a declared-low-MOI screen now compares
  each perturbation against the non-targeting cells, the contrast SCEPTRE documents
  for low-MOI data and the one PerTurbo was already using. **SCEPTRE p-values and
  effect sizes on low-MOI screens will differ from every run made before this
  change.** High-MOI screens are unaffected, and PerTurbo's behaviour is unchanged
  in both cases. To reproduce the old low-MOI SCEPTRE behaviour set
  `INFERENCE_control_group = 'complement'` (which also switches PerTurbo to the
  all-cells pool) or override SCEPTRE alone with
  `INFERENCE_SCEPTRE_control_group`.
- Asking for `nt_cells` on a high-MOI analysis now fails with an explicit message
  instead of silently substituting `complement`.

### `Deprecated`

- `INFERENCE_PERTURBO_CRT_POOL` and `INFERENCE_SCEPTRE_control_group` remain for
  back-compatibility as per-method overrides. Their historical defaults
  (`from-moi`, `complement`) count as unset; any other value wins for that method
  alone and the pipeline logs that the two methods are deliberately inconsistent.

## v1.0.0 - 2025-AUG-19

- Initial release of IGVF Perturb-seq Pipeline
- Core pipeline functionality
- Documentation and examples

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
