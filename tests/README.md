# CRISPR Pipeline Tests

This directory contains tests for the sceptre functionality updates made to the CRISPR Pipeline.

## Directory Structure

```
tests/
├── README.md                    # This file
├── run_tests.sh                 # Test runner script
├── data/                        # Test datasets from sceptreIGVF
│   ├── mudata_guide_assignment_gasperini.rda
│   ├── mudata_guide_assignment_papalexi.rda
│   ├── mudata_inference_gasperini.rda
│   └── mudata_inference_papalexi.rda
├── unit/                        # Unit tests for individual functions
│   ├── test_grna_assignment.R   # Tests for gRNA assignment with new parameters
│   └── test_inference.R         # Tests for updated inference functionality
└── integration/                 # Integration tests
    └── test_workflow.R          # End-to-end workflow tests
```

## What's Being Tested

### gRNA Assignment Tests (`unit/test_grna_assignment.R`)
- ✅ Default parameter behavior (backward compatibility)
- ✅ Custom `probability_threshold` parameter
- ✅ Custom `n_em_rep` parameter  
- ✅ Both custom parameters together
- ✅ Parameter type conversion (string to numeric/integer)

### Inference Tests (`unit/test_inference.R`)
- ✅ `convert_mudata_to_sceptre_object_v1()` function
- ✅ Basic inference with default parameters
- ✅ Inference with additional sceptre parameters
- ✅ Sceptre object validity
- ✅ Complete inference workflow
- ✅ Result consistency across runs

### Integration Tests (`integration/test_workflow.R`)
- ✅ Complete workflow: assignment → inference (default parameters)
- ✅ Complete workflow with custom parameters
- ✅ Data integrity through workflow
- ✅ Error handling for edge cases
- ✅ Performance validation

## Requirements

The following R packages are required to run the tests:
- `testthat` - Testing framework
- `MuData` - Multi-omics data handling
- `MultiAssayExperiment` - Multi-assay data structures
- `sceptre` - Single-cell perturbation analysis
- `Matrix` - Sparse matrix operations
- `dplyr` - Data manipulation

## Running Tests

### Run All Tests
```bash
# From the tests directory
./run_tests.sh
```

### Run Individual Tests
```bash
# From the tests directory
Rscript unit/test_grna_assignment.R
Rscript unit/test_inference.R
Rscript integration/test_workflow.R
```

## Test Data

The test data is copied from the sceptreIGVF R package and includes:
- **Gasperini dataset**: CRISPR screen data for testing
- **Papalexi dataset**: Alternative CRISPR screen data

Both datasets include:
- Gene expression data
- Guide RNA data  
- Metadata for analysis parameters

## Expected Behavior

All tests should pass if the sceptre parameter updates are working correctly:

- **Default parameters**: Should maintain backward compatibility
- **Custom parameters**: Should properly pass parameters to underlying sceptre functions
- **Data integrity**: Original data should be preserved through the workflow
- **Results**: Should produce valid statistical results (p-values, log2 fold changes)

## Troubleshooting

If tests fail:
1. Check that all required R packages are installed
2. Verify that sceptre package is properly installed and functional
3. Check that test data files are present in `tests/data/`
4. Review error messages for specific function failures

The test framework uses `testthat` which provides detailed error reporting for debugging.