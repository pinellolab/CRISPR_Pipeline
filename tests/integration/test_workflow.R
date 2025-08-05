#!/usr/bin/env Rscript

# Integration test for complete gRNA assignment + inference workflow
# Tests the full pipeline: assignment -> inference

# Load required libraries
library(testthat)
library(MuData)
library(MultiAssayExperiment)

# Source both functions (set flag to prevent command execution)
.sourced_from_test <- TRUE
source("../../bin/assign_grnas_sceptre.R")
source("../../bin/inference_sceptre.R")

cat("Starting integration workflow tests...\n")

# Load test data (use assignment data since we'll run assignment first)
load("../data/mudata_guide_assignment_gasperini.rda")

# Integration Test 1: Default parameters workflow
test_that("complete workflow works with default parameters", {
  
  cat("Testing complete workflow with default parameters...\n")
  
  # Step 1: gRNA assignment
  assignment_result <- assign_grnas_sceptre_v1(
    mudata_guide_assignment_gasperini,
    probability_threshold = "default",
    n_em_rep = "default"
  )
  
  expect_true("mudata" %in% names(assignment_result))
  expect_true("guide_assignment" %in% names(assignment_result))
  
  # Get the mudata with guide assignments
  mudata_with_assignments <- assignment_result$mudata
  
  # Verify guide assignment was added
  expect_true("guide_assignment" %in% SummarizedExperiment::assayNames(mudata_with_assignments[['guide']]))
  
  # Step 2: Create some test pairs for inference
  # (In real pipeline this would come from elsewhere)
  n_genes <- nrow(mudata_with_assignments[['gene']])
  
  # Get actual guide targets from the data
  guide_targets <- unique(SingleCellExperiment::rowData(mudata_with_assignments[['guide']])$intended_target_name)
  # Use first non "non-targeting" target
  test_target <- guide_targets[guide_targets != "non-targeting"][1]
  
  # Create a small set of test pairs
  test_pairs <- data.frame(
    gene_id = rownames(mudata_with_assignments[['gene']])[1:min(10, n_genes)],
    intended_target_name = rep(test_target, min(10, n_genes))
  )
  
  MultiAssayExperiment::metadata(mudata_with_assignments)$pairs_to_test <- test_pairs
  
  # Step 3: Run inference
  inference_result <- inference_sceptre_m(mudata_with_assignments)
  
  expect_true("mudata" %in% names(inference_result))
  expect_true("test_results" %in% names(inference_result))
  
  # Validate inference results
  test_results <- inference_result$test_results
  expect_true(is.data.frame(test_results))
  expect_true("p_value" %in% colnames(test_results))
  expect_true("log2_fc" %in% colnames(test_results))
  
  cat("✓ Default parameters workflow test passed\n")
})

# Integration Test 2: Custom parameters workflow
test_that("complete workflow works with custom parameters", {
  
  cat("Testing complete workflow with custom parameters...\n")
  
  # Step 1: gRNA assignment with custom parameters
  assignment_result <- assign_grnas_sceptre_v1(
    mudata_guide_assignment_gasperini,
    probability_threshold = "0.8",
    n_em_rep = "10"
  )
  
  expect_true("guide_assignment" %in% names(assignment_result))
  
  mudata_with_assignments <- assignment_result$mudata
  
  # Step 2: Create test pairs
  guide_targets <- unique(SingleCellExperiment::rowData(mudata_with_assignments[['guide']])$intended_target_name)
  test_target <- guide_targets[guide_targets != "non-targeting"][1]
  
  test_pairs <- data.frame(
    gene_id = rownames(mudata_with_assignments[['gene']])[1:5],
    intended_target_name = rep(test_target, 5)
  )
  
  MultiAssayExperiment::metadata(mudata_with_assignments)$pairs_to_test <- test_pairs
  
  # Step 3: Run inference with custom parameters
  inference_result <- inference_sceptre_m(
    mudata_with_assignments,
    side = "both",
    grna_integration_strategy = "union"
  )
  
  expect_true("test_results" %in% names(inference_result))
  
  # Check that we got reasonable results
  test_results <- inference_result$test_results
  expect_true(nrow(test_results) == nrow(test_pairs))
  expect_true(all(test_results$p_value >= 0 & test_results$p_value <= 1, na.rm = TRUE))
  
  cat("✓ Custom parameters workflow test passed\n")
})

# Integration Test 3: Data consistency through workflow
test_that("data integrity is maintained through workflow", {
  
  cat("Testing data integrity through workflow...\n")
  
  original_mudata <- mudata_guide_assignment_gasperini
  
  # Run assignment
  assignment_result <- assign_grnas_sceptre_v1(original_mudata)
  processed_mudata <- assignment_result$mudata
  
  # Check that original data is preserved
  expect_equal(dim(processed_mudata[['gene']]), dim(original_mudata[['gene']]))
  expect_equal(dim(processed_mudata[['guide']]), dim(original_mudata[['guide']]))
  
  # Check that guide assignment was added
  original_assays <- SummarizedExperiment::assayNames(original_mudata[['guide']])
  processed_assays <- SummarizedExperiment::assayNames(processed_mudata[['guide']])
  
  # Original assays might be NULL, so handle that case
  original_count <- ifelse(is.null(original_assays), 0, length(original_assays))
  processed_count <- length(processed_assays)
  
  # Should have at least one more assay (guide_assignment)
  expect_true(processed_count >= original_count + 1)
  expect_true("guide_assignment" %in% processed_assays)
  
  cat("✓ Data integrity test passed\n")
})

# Integration Test 4: Error handling
test_that("workflow handles edge cases gracefully", {
  
  cat("Testing error handling...\n")
  
  # Test with minimal test pairs
  assignment_result <- assign_grnas_sceptre_v1(mudata_guide_assignment_gasperini)
  mudata_with_assignments <- assignment_result$mudata
  
  # Create minimal test pairs (just 1 pair)
  guide_targets <- unique(SingleCellExperiment::rowData(mudata_with_assignments[['guide']])$intended_target_name)
  test_target <- guide_targets[guide_targets != "non-targeting"][1]
  
  minimal_pairs <- data.frame(
    gene_id = rownames(mudata_with_assignments[['gene']])[1],
    intended_target_name = test_target
  )
  
  MultiAssayExperiment::metadata(mudata_with_assignments)$pairs_to_test <- minimal_pairs
  
  # This should still work
  expect_no_error({
    inference_result <- inference_sceptre_m(mudata_with_assignments)
  })
  
  cat("✓ Error handling test passed\n")
})

# Integration Test 5: Performance check
test_that("workflow completes in reasonable time", {
  
  cat("Testing workflow performance...\n")
  
  start_time <- Sys.time()
  
  # Run a minimal workflow
  assignment_result <- assign_grnas_sceptre_v1(mudata_guide_assignment_gasperini)
  mudata_with_assignments <- assignment_result$mudata
  
  # Small test set
  guide_targets <- unique(SingleCellExperiment::rowData(mudata_with_assignments[['guide']])$intended_target_name)
  test_target <- guide_targets[guide_targets != "non-targeting"][1]
  
  test_pairs <- data.frame(
    gene_id = rownames(mudata_with_assignments[['gene']])[1:3],
    intended_target_name = rep(test_target, 3)
  )
  
  MultiAssayExperiment::metadata(mudata_with_assignments)$pairs_to_test <- test_pairs
  
  inference_result <- inference_sceptre_m(mudata_with_assignments)
  
  end_time <- Sys.time()
  runtime <- as.numeric(end_time - start_time, units = "secs")
  
  cat(sprintf("Workflow completed in %.2f seconds\n", runtime))
  
  # Should complete within reasonable time (adjust threshold as needed)
  expect_true(runtime < 300) # 5 minutes max for test
  
  cat("✓ Performance test passed\n")
})

cat("\nAll integration tests completed successfully! ✓\n")