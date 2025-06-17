#!/bin/bash

# Test helper functions for enhanced integration testing

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test tracking
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Initialize test environment
init_test_env() {
    export COVTOBED_QUIET=1
    TESTS_RUN=0
    TESTS_PASSED=0
    TESTS_FAILED=0
}

# Print colored output
print_status() {
    local status="$1"
    local message="$2"
    
    case "$status" in
        "PASS")
            echo -e "${GREEN}[PASS]${NC} $message"
            ;;
        "FAIL")
            echo -e "${RED}[FAIL]${NC} $message"
            ;;
        "INFO")
            echo -e "${BLUE}[INFO]${NC} $message"
            ;;
        "WARN")
            echo -e "${YELLOW}[WARN]${NC} $message"
            ;;
    esac
}

# Run a test with enhanced error reporting
run_test() {
    local test_name="$1"
    local test_command="$2"
    local expected_exit_code="${3:-0}"
    
    TESTS_RUN=$((TESTS_RUN + 1))
    
    echo -n "Testing $test_name: "
    
    # Create temporary files for output capture
    local stdout_file=$(mktemp)
    local stderr_file=$(mktemp)
    
    # Run the test command
    eval "$test_command" > "$stdout_file" 2> "$stderr_file"
    local actual_exit_code=$?
    
    # Check exit code
    if [ "$actual_exit_code" -eq "$expected_exit_code" ]; then
        print_status "PASS" "$test_name"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        local result=0
    else
        print_status "FAIL" "$test_name (exit code: expected $expected_exit_code, got $actual_exit_code)"
        echo "STDOUT:"
        cat "$stdout_file" | sed 's/^/  /'
        echo "STDERR:"
        cat "$stderr_file" | sed 's/^/  /'
        TESTS_FAILED=$((TESTS_FAILED + 1))
        local result=1
    fi
    
    # Cleanup
    rm -f "$stdout_file" "$stderr_file"
    return $result
}

# Assert functions
assert_equals() {
    local expected="$1"
    local actual="$2"
    local message="${3:-Assertion failed}"
    
    if [ "$expected" = "$actual" ]; then
        return 0
    else
        print_status "FAIL" "$message (expected: '$expected', actual: '$actual')"
        return 1
    fi
}

assert_not_equals() {
    local not_expected="$1"
    local actual="$2"
    local message="${3:-Assertion failed}"
    
    if [ "$not_expected" != "$actual" ]; then
        return 0
    else
        print_status "FAIL" "$message (should not equal: '$not_expected')"
        return 1
    fi
}

assert_not_empty() {
    local value="$1"
    local message="${2:-Value should not be empty}"
    
    if [ -n "$value" ]; then
        return 0
    else
        print_status "FAIL" "$message"
        return 1
    fi
}

assert_file_exists() {
    local file="$1"
    local message="${2:-File should exist}"
    
    if [ -f "$file" ]; then
        return 0
    else
        print_status "FAIL" "$message: $file"
        return 1
    fi
}

assert_line_count() {
    local expected_lines="$1"
    local output="$2"
    local message="${3:-Line count assertion failed}"
    
    local actual_lines=$(echo "$output" | wc -l)
    assert_equals "$expected_lines" "$actual_lines" "$message"
}

assert_contains() {
    local haystack="$1"
    local needle="$2"
    local message="${3:-String should contain pattern}"
    
    if echo "$haystack" | grep -q "$needle"; then
        return 0
    else
        print_status "FAIL" "$message (looking for '$needle')"
        return 1
    fi
}

# Performance testing helpers
measure_time() {
    local command="$1"
    local iterations="${2:-1}"
    
    local total_time=0
    for ((i=1; i<=iterations; i++)); do
        local start_time=$(date +%s.%N)
        eval "$command" > /dev/null 2>&1
        local end_time=$(date +%s.%N)
        local duration=$(echo "$end_time - $start_time" | bc -l)
        total_time=$(echo "$total_time + $duration" | bc -l)
    done
    
    local avg_time=$(echo "scale=3; $total_time / $iterations" | bc -l)
    echo "$avg_time"
}

# Memory usage measurement
measure_memory() {
    local command="$1"
    local max_memory=0
    
    # Run command in background and monitor memory
    eval "$command" &
    local pid=$!
    
    while kill -0 "$pid" 2>/dev/null; do
        if [ -f "/proc/$pid/status" ]; then
            local current_memory=$(grep VmRSS /proc/$pid/status | awk '{print $2}')
            if [ "$current_memory" -gt "$max_memory" ]; then
                max_memory=$current_memory
            fi
        fi
        sleep 0.1
    done
    
    wait "$pid"
    echo "$max_memory"
}

# Summary reporting
print_test_summary() {
    echo ""
    echo "=========================================="
    echo "Test Summary:"
    echo "  Total tests: $TESTS_RUN"
    print_status "PASS" "Passed: $TESTS_PASSED"
    if [ "$TESTS_FAILED" -gt 0 ]; then
        print_status "FAIL" "Failed: $TESTS_FAILED"
        return 1
    else
        print_status "PASS" "All tests passed!"
        return 0
    fi
}