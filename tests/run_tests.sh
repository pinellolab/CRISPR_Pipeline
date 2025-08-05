#!/bin/bash

# Test runner script for CRISPR Pipeline sceptre functionality tests

echo "=========================================="
echo "Running CRISPR Pipeline Tests"
echo "=========================================="

# Check if required libraries are available
echo "Checking R environment..."
Rscript -e "
required_packages <- c('testthat', 'MuData', 'MultiAssayExperiment', 'sceptre', 'Matrix', 'dplyr')
missing_packages <- required_packages[!sapply(required_packages, require, character.only = TRUE, quietly = TRUE)]
if(length(missing_packages) > 0) {
  cat('Missing required packages:', paste(missing_packages, collapse = ', '), '\n')
  cat('Please install missing packages before running tests.\n')
  quit(status = 1)
} else {
  cat('All required packages are available.\n')
}
"

if [ $? -ne 0 ]; then
    echo "❌ Missing required R packages. Please install them first."
    exit 1
fi

# Set up test environment
cd "$(dirname "$0")"
TEST_DIR=$(pwd)

echo ""
echo "Test directory: $TEST_DIR"
echo ""

# Initialize test results
PASSED=0
FAILED=0

# Function to run a test and track results
run_test() {
    local test_file=$1
    local test_name=$2
    
    echo "----------------------------------------"
    echo "Running: $test_name"
    echo "----------------------------------------"
    
    # Change to the directory containing the test file to ensure relative paths work
    local test_dir=$(dirname "$test_file")
    local test_filename=$(basename "$test_file")
    
    if (cd "$test_dir" && Rscript "$test_filename"); then
        echo "✅ $test_name PASSED"
        ((PASSED++))
    else
        echo "❌ $test_name FAILED"
        ((FAILED++))
    fi
    echo ""
}

# Run unit tests
echo "Running Unit Tests..."
echo "===================="

run_test "unit/test_grna_assignment.R" "gRNA Assignment Tests"
run_test "unit/test_inference.R" "Inference Tests"

# Run integration tests
echo "Running Integration Tests..."
echo "============================"

run_test "integration/test_workflow.R" "Workflow Integration Tests"

# Summary
echo "=========================================="
echo "TEST SUMMARY"
echo "=========================================="
echo "Passed: $PASSED"
echo "Failed: $FAILED"
echo "Total:  $((PASSED + FAILED))"

if [ $FAILED -eq 0 ]; then
    echo ""
    echo "🎉 All tests passed!"
    exit 0
else
    echo ""
    echo "💥 Some tests failed. Please review the output above."
    exit 1
fi