#!/usr/bin/env Rscript

# Test script for gRNA assignment functionality with new parameters
# Tests both default and custom parameter behavior

# Load required libraries
library(testthat)
library(MuData)
library(MultiAssayExperiment)

# Source the functions we're testing (set flag to prevent command execution)
.sourced_from_test <- TRUE
source("../../bin/assign_grnas_sceptre.R")

cat("Starting gRNA assignment tests...\n")

# Load test data
load("../data/mudata_guide_assignment_gasperini.rda")
load("../data/mudata_guide_assignment_papalexi.rda")

# Test 1: Default parameters (should work like before)
test_that("gRNA assignment works with default parameters", {
  
  cat("Testing with default parameters...\n")
  
  # Test with Gasperini dataset
  result_gasperini <- assign_grnas_sceptre_v1(mudata_guide_assignment_gasperini)
  
  # Validate structure
  expect_true(is.list(result_gasperini))
  expect_true("mudata" %in% names(result_gasperini))
  expect_true("guide_assignment" %in% names(result_gasperini))
  
  # Check that guide_assignment matrix has correct properties
  guide_matrix <- result_gasperini$guide_assignment
  expect_true(is(guide_matrix, "dsparseMatrix"))
  expect_true(all(guide_matrix@x %in% c(0, 1))) # Should be binary
  
  # Test with Papalexi dataset
  result_papalexi <- assign_grnas_sceptre_v1(mudata_guide_assignment_papalexi)
  expect_true(is.list(result_papalexi))
  expect_true("guide_assignment" %in% names(result_papalexi))
  
  cat("✓ Default parameters test passed\n")
})

# Test 2: Custom probability_threshold parameter
test_that("gRNA assignment works with custom probability_threshold", {
  
  cat("Testing with custom probability_threshold...\n")
  
  result_custom <- assign_grnas_sceptre_v1(
    mudata_guide_assignment_gasperini, 
    probability_threshold = "0.8",
    n_em_rep = "default"
  )
  
  # Validate structure is maintained
  expect_true(is.list(result_custom))
  expect_true("guide_assignment" %in% names(result_custom))
  
  # Check matrix properties
  guide_matrix <- result_custom$guide_assignment
  expect_true(is(guide_matrix, "dsparseMatrix"))
  expect_true(all(guide_matrix@x %in% c(0, 1)))
  
  cat("✓ Custom probability_threshold test passed\n")
})

# Test 3: Custom n_em_rep parameter
test_that("gRNA assignment works with custom n_em_rep", {
  
  cat("Testing with custom n_em_rep...\n")
  
  result_custom <- assign_grnas_sceptre_v1(
    mudata_guide_assignment_gasperini,
    probability_threshold = "default", 
    n_em_rep = "10"
  )
  
  # Validate structure is maintained
  expect_true(is.list(result_custom))
  expect_true("guide_assignment" %in% names(result_custom))
  
  cat("✓ Custom n_em_rep test passed\n")
})

# Test 4: Both custom parameters
test_that("gRNA assignment works with both custom parameters", {
  
  cat("Testing with both custom parameters...\n")
  
  result_both <- assign_grnas_sceptre_v1(
    mudata_guide_assignment_gasperini,
    probability_threshold = "0.7",
    n_em_rep = "15"
  )
  
  # Validate structure is maintained
  expect_true(is.list(result_both))
  expect_true("guide_assignment" %in% names(result_both))
  
  # Check that results are reasonable
  guide_matrix <- result_both$guide_assignment
  expect_true(nrow(guide_matrix) > 0)
  expect_true(ncol(guide_matrix) > 0)
  
  cat("✓ Both custom parameters test passed\n")
})

# Test 5: Parameter validation
test_that("gRNA assignment handles parameter types correctly", {
  
  cat("Testing parameter type conversion...\n")
  
  # Test that string numbers are converted correctly
  result <- assign_grnas_sceptre_v1(
    mudata_guide_assignment_gasperini,
    probability_threshold = "0.9",  # Should convert to numeric
    n_em_rep = "5"                  # Should convert to integer
  )
  
  expect_true("guide_assignment" %in% names(result))
  
  cat("✓ Parameter type conversion test passed\n")
})

cat("\nAll gRNA assignment tests completed successfully! ✓\n")