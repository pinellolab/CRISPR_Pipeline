# Pull Request

## Add sceptre parameter exposure and fix inference multicollinearity

### Summary

This PR enhances the sceptre integration in the Nextflow pipeline by:
- Exposing advanced gRNA assignment parameters (`probability_threshold` and `n_em_rep`)
- Switching to sceptre's default formula object for inference, while excluding gRNA technical covariates to avoid multicollinearity issues.
- Switching inference to use `convert_mudata_to_sceptre_object_v1()`, for consistency with the guide assignment functionality.

### Changes Made

#### gRNA Assignment Enhancements
- **Function Updates**: Modified `assign_grnas_sceptre_v1()` to accept `probability_threshold` and `n_em_rep` parameters
- **Parameter Handling**: Implemented "default" fallback behavior - when set to "default", parameters are omitted to inherit sceptre's defaults
- **Pipeline Integration**: Updated `guide_assignment_sceptre` Nextflow process to pass parameters from pipeline configuration
- **Configuration**: Added new parameters to `nextflow.config`:
  - `GUIDE_ASSIGNMENT_SCEPTRE_probability_threshold = 'default'`
  - `GUIDE_ASSIGNMENT_SCEPTRE_n_em_rep = 'default'`

#### Inference Function Improvements
- **Local Function Usage**: Switched from `sceptreIGVF::convert_mudata_to_sceptre_object()` to local `convert_mudata_to_sceptre_object_v1()`
- **Multicollinearity Fix**: Added custom formula construction that excludes gRNA technical covariates (`grna_n_nonzero` and `grna_n_umis`) which become identical when using binary gRNA assignments
- **Formula Handling**: Uses `sceptre:::auto_construct_formula_object()` with `include_grna_covariates = FALSE` to prevent redundant covariate errors
- **Parameter Cleanup**: Removed `INFERENCE_SCEPTRE_formula_object` from pipeline configuration
- **Process Simplification**: Updated `inference_sceptre` Nextflow process to remove unnecessary parameter passing

### Files Changed

- `bin/assign_grnas_sceptre.R` - Added parameter support to assignment function
- `bin/inference_sceptre.R` - Fixed multicollinearity issue and added local conversion
- `modules/local/guide_assignment_sceptre/main.nf` - Updated to pass new parameters
- `modules/local/inference_sceptre/main.nf` - Simplified parameter structure
- `nextflow.config` - Added new parameters and removed deprecated ones
