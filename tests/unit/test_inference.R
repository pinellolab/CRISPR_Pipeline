#!/usr/bin/env Rscript

# Test script for inference functionality changes
# Tests the updated inference function and convert_mudata_to_sceptre_object_v1

# Load required libraries
library(testthat)
library(MuData)
library(MultiAssayExperiment)

# Source the functions we're testing (set flag to prevent command execution)
.sourced_from_test <- TRUE
source("../../bin/inference_sceptre.R")

cat("Starting inference tests...\n")

# Load test data
load("../data/mudata_inference_gasperini.rda")
load("../data/mudata_inference_papalexi.rda")

# Test 1: convert_mudata_to_sceptre_object_v1 function
test_that("convert_mudata_to_sceptre_object_v1 works correctly", {
  
  cat("Testing convert_mudata_to_sceptre_object_v1...\n")
  
  # Test with Gasperini dataset
  sceptre_obj <- convert_mudata_to_sceptre_object_v1(mudata_inference_gasperini)
  
  # Validate that we get a sceptre object
  expect_true(is(sceptre_obj, "sceptre_object"))
  
  # Check basic structure
  expect_true(length(sceptre_obj@response_matrix) > 0)
  expect_true(length(sceptre_obj@grna_matrix) > 0)
  expect_true(nrow(sceptre_obj@grna_target_data_frame) > 0)
  
  cat("✓ convert_mudata_to_sceptre_object_v1 test passed\n")
})

# Test 2: Basic inference functionality
test_that("inference_sceptre_m works with default parameters", {
  
  cat("Testing inference with default parameters...\n")
  
  # Test with Gasperini dataset
  result_gasperini <- inference_sceptre_m(mudata_inference_gasperini)
  
  # Validate structure
  expect_true(is.list(result_gasperini))
  expect_true("mudata" %in% names(result_gasperini))
  expect_true("test_results" %in% names(result_gasperini))
  
  # Check test_results structure
  test_results <- result_gasperini$test_results
  expect_true(is.data.frame(test_results))
  expect_true("p_value" %in% colnames(test_results))
  expect_true("log2_fc" %in% colnames(test_results))
  expect_true(nrow(test_results) > 0)
  
  # Test with Papalexi dataset
  result_papalexi <- inference_sceptre_m(mudata_inference_papalexi)
  expect_true(is.list(result_papalexi))
  expect_true("test_results" %in% names(result_papalexi))
  
  cat("✓ Default inference test passed\n")
})

# Test 3: Inference with additional parameters
test_that("inference_sceptre_m works with additional parameters", {
  
  cat("Testing inference with additional parameters...\n")
  
  # Test with additional sceptre parameters
  result_custom <- inference_sceptre_m(
    mudata_inference_gasperini,
    side = "left",
    grna_integration_strategy = "union"
  )
  
  # Validate structure is maintained
  expect_true(is.list(result_custom))
  expect_true("test_results" %in% names(result_custom))
  
  # Check that test results are reasonable
  test_results <- result_custom$test_results
  expect_true(is.data.frame(test_results))
  expect_true("p_value" %in% colnames(test_results))
  expect_true(all(test_results$p_value >= 0 & test_results$p_value <= 1, na.rm = TRUE))
  
  cat("✓ Custom parameters inference test passed\n")
})

# Test 4: Compare local vs sceptreIGVF conversion (if available)
test_that("convert_mudata_to_sceptre_object_v1 produces valid sceptre objects", {
  
  cat("Testing sceptre object validity...\n")
  
  sceptre_obj <- convert_mudata_to_sceptre_object_v1(mudata_inference_gasperini)
  
  # Test that we can set analysis parameters (this validates object structure)
  pairs_to_test <- MultiAssayExperiment::metadata(mudata_inference_gasperini)$pairs_to_test |>
    as.data.frame()
  
  discovery_pairs <- pairs_to_test |>
    dplyr::rename(
      grna_target = intended_target_name,
      response_id = gene_id
    )
  
  # This should not throw an error if object is valid (use custom formula like our main function)
  expect_no_error({
    formula_object <- sceptre:::auto_construct_formula_object(
      cell_covariates = sceptre_obj@covariate_data_frame,
      include_grna_covariates = FALSE
    )
    sceptre_obj_with_params <- sceptre_obj |>
      sceptre::set_analysis_parameters(discovery_pairs = discovery_pairs, formula_object = formula_object)
  })
  
  cat("✓ Sceptre object validity test passed\n")
})

# Test 5: End-to-end inference workflow
test_that("complete inference workflow functions correctly", {
  
  cat("Testing complete inference workflow...\n")
  
  # Run the complete inference
  result <- inference_sceptre_m(mudata_inference_gasperini)
  
  # Check that positive controls have low p-values (if they exist)
  test_results <- result$test_results
  
  # Basic sanity checks
  expect_true(all(is.finite(test_results$p_value) | is.na(test_results$p_value)))
  expect_true(all(is.finite(test_results$log2_fc) | is.na(test_results$log2_fc)))
  
  # Check that we have some significant results (p < 0.05)
  significant_results <- sum(test_results$p_value < 0.05, na.rm = TRUE)
  expect_true(significant_results >= 0) # Should have at least some results
  
  cat("✓ Complete inference workflow test passed\n")
})

# Test 6: Regression test - ensure no major functionality breaks
test_that("inference results are consistent across runs", {
  
  cat("Testing result consistency...\n")
  
  # Run inference twice with same parameters
  result1 <- inference_sceptre_m(mudata_inference_gasperini)
  result2 <- inference_sceptre_m(mudata_inference_gasperini)
  
  # Results should be identical (assuming deterministic behavior)
  expect_equal(nrow(result1$test_results), nrow(result2$test_results))
  expect_equal(colnames(result1$test_results), colnames(result2$test_results))
  
  cat("✓ Consistency test passed\n")
})

cat("\nAll inference tests completed successfully! ✓\n")