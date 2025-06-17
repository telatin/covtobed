#!/bin/bash

# Enhanced integration tests for covtobed
# This script provides more comprehensive testing than the original test.sh

set -eou pipefail

# Get script directory and source helpers
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/test_helpers.sh"

# Change to project root
cd "$SCRIPT_DIR/../.."

# Initialize test environment
init_test_env

print_status "INFO" "Starting enhanced integration tests for covtobed"

# Check if binary exists
if [ ! -e "./covtobed" ]; then
    print_status "FAIL" "covtobed binary not found. Please compile first."
    exit 1
fi

# Test 1: Basic version check
run_test "version_check" \
    "./covtobed --version | head -1" \
    0

# Test 2: Help output
run_test "help_output" \
    "./covtobed --help" \
    0

# Test 3: Basic BAM processing
run_test "basic_bam_processing" \
    "./covtobed test/demo.bam | wc -l | grep -q 202" \
    0

# Test 4: Sorted vs unsorted BAM
run_test "sorted_bam_accepted" \
    "./covtobed test/mini-sorted.bam" \
    0

run_test "unsorted_bam_rejected" \
    "./covtobed test/mini-unsorted.bam" \
    1

# Test 5: Coverage filtering
run_test "minimum_coverage_filter" \
    "./covtobed -m 15 test/demo.bam | wc -l | grep -q 12" \
    0

run_test "maximum_coverage_filter" \
    "./covtobed -x 5 test/demo.bam | head -10" \
    0

# Test 6: Physical coverage
run_test "physical_coverage" \
    "./covtobed --physical-coverage test/mp.bam | wc -l | grep -q 136" \
    0

# Test 7: Stranded output
run_test "stranded_output" \
    "./covtobed --output-strands test/demo.bam | cut -f 5 | sort -u | wc -l | grep -q 10" \
    0

# Test 8: Different output formats
run_test "bed_format_output" \
    "./covtobed --format bed test/demo.bam | head -1" \
    0

run_test "counts_format_output" \
    "./covtobed --format counts test/demo.bam | grep '>' | wc -l | grep -q 2" \
    0

run_test "counts_format_lines" \
    "./covtobed --format counts test/demo.bam | grep -v '>' | wc -l | grep -q 202" \
    0

# Test 9: Length filtering
run_test "minimum_length_filter" \
    "./covtobed -l 50 test/demo.bam" \
    0

# Test 10: Quality filtering
run_test "mapping_quality_filter" \
    "./covtobed -q 20 test/demo.bam" \
    0

# Test 11: Invalid alignment filtering
run_test "discard_invalid_alignments" \
    "./covtobed -d test/filtered.bam | wc -l | grep -q 2" \
    0

run_test "keep_invalid_alignments" \
    "./covtobed test/filtered.bam | wc -l | grep -q 6" \
    0

# Test 12: Output consistency check
print_status "INFO" "Running output consistency checks..."

TEMP_OUTPUT=$(mktemp)
./covtobed test/demo.bam > "$TEMP_OUTPUT"

if diff -q "$TEMP_OUTPUT" test/output.bed > /dev/null; then
    print_status "PASS" "Output matches reference"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    print_status "FAIL" "Output differs from reference"
    echo "Differences:"
    diff "$TEMP_OUTPUT" test/output.bed | head -10
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

rm -f "$TEMP_OUTPUT"
TESTS_RUN=$((TESTS_RUN + 1))

# Test 13: Synthetic BAM validation
TEMP_MOCK_OUTPUT=$(mktemp)
./covtobed -m 1 test/mock.bam > "$TEMP_MOCK_OUTPUT"

if diff -q "$TEMP_MOCK_OUTPUT" test/mock.bed > /dev/null; then
    print_status "PASS" "Synthetic BAM output correct"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    print_status "FAIL" "Synthetic BAM output incorrect"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi

rm -f "$TEMP_MOCK_OUTPUT"
TESTS_RUN=$((TESTS_RUN + 1))

# Test 14: Edge cases
print_status "INFO" "Testing edge cases..."

# Empty BAM handling (if available)
if [ -f "test/empty.bam" ]; then
    run_test "empty_bam_handling" \
        "./covtobed test/empty.bam" \
        0
fi

# Very small BAM
run_test "minimal_bam_processing" \
    "./covtobed test/nano.bam" \
    0

# Test 15: Performance test (basic)
print_status "INFO" "Running basic performance test..."

PERF_TIME=$(measure_time "./covtobed test/demo.bam > /dev/null" 3)
print_status "INFO" "Average processing time: ${PERF_TIME}s"

if (( $(echo "$PERF_TIME < 1.0" | bc -l) )); then
    print_status "PASS" "Performance acceptable (< 1.0s)"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    print_status "WARN" "Performance slower than expected (${PERF_TIME}s)"
fi
TESTS_RUN=$((TESTS_RUN + 1))

# Test 16: Memory usage test (basic)
print_status "INFO" "Running basic memory test..."

if command -v valgrind >/dev/null 2>&1; then
    VALGRIND_OUTPUT=$(mktemp)
    if valgrind --tool=massif --stacks=yes --massif-out-file="$VALGRIND_OUTPUT" \
       ./covtobed test/demo.bam > /dev/null 2>&1; then
        print_status "PASS" "Memory test completed (check $VALGRIND_OUTPUT for details)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        print_status "FAIL" "Memory test failed"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
    TESTS_RUN=$((TESTS_RUN + 1))
else
    print_status "INFO" "Valgrind not available, skipping memory test"
fi

# Test 17: Strand-specific edge cases
run_test "strand_adjacent_intervals" \
    "[ \$(./covtobed test/stranded.bam | wc -l) -eq 1 ] && [ \$(./covtobed --output-strands test/stranded.bam | wc -l) -eq 2 ]" \
    0

# Test 18: Artificial coverage validation
print_status "INFO" "Validating artificial coverage patterns..."

ARTIFICIAL_TEST_PASSED=true
./covtobed test/test_cov.bam -m 1 | cut -f 1,4 | while read LINE; do
    if ! echo "$LINE" | perl -ne '($exp, $cov)=split /\s+/, $_; 
        if ("$exp" ne "${cov}X") {
            print STDERR "FAIL: $exp != ${cov}X\n"; 
            exit 1;
        }'; then
        ARTIFICIAL_TEST_PASSED=false
        break
    fi
done

if [ "$ARTIFICIAL_TEST_PASSED" = true ]; then
    print_status "PASS" "Artificial coverage validation"
    TESTS_PASSED=$((TESTS_PASSED + 1))
else
    print_status "FAIL" "Artificial coverage validation"
    TESTS_FAILED=$((TESTS_FAILED + 1))
fi
TESTS_RUN=$((TESTS_RUN + 1))

# Print summary
print_test_summary