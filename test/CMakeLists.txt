cmake_minimum_required(VERSION 3.10)
project(covtobed_tests)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add subdirectories
add_subdirectory(unit)
add_subdirectory(benchmark)
add_subdirectory(utils)

# Create a target that runs all tests
add_custom_target(run_all_tests
    COMMAND echo "Running unit tests..."
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/unit/unit_tests
    COMMAND echo "Running integration tests..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_enhanced.sh
    DEPENDS unit_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
)

# Target for running benchmarks
add_custom_target(run_all_benchmarks
    COMMAND echo "Running benchmarks..."
    COMMAND ${CMAKE_CURRENT_BINARY_DIR}/benchmark/benchmark_coverage
    DEPENDS benchmark_coverage
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/benchmark
)

# Generate synthetic test data
add_custom_target(generate_test_data
    COMMAND echo "Generating synthetic test data..."
    COMMAND mkdir -p ${CMAKE_CURRENT_SOURCE_DIR}/data/synthetic
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_test_generator.sh ${CMAKE_CURRENT_SOURCE_DIR}/data/synthetic/coverage_pattern.sam coverage 100 200 10 || echo "Warning: Failed to generate coverage_pattern.sam"
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_test_generator.sh ${CMAKE_CURRENT_SOURCE_DIR}/data/synthetic/stranded_test.sam stranded || echo "Warning: Failed to generate stranded_test.sam"
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/utils/simple_test_generator.sh ${CMAKE_CURRENT_SOURCE_DIR}/data/synthetic/single_read.sam single || echo "Warning: Failed to generate single_read.sam"
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# Comprehensive test target
add_custom_target(test_comprehensive
    COMMAND echo "=== Running comprehensive test suite ==="
    COMMAND echo "1. Generating test data..."
    COMMAND $(MAKE) generate_test_data
    COMMAND echo "2. Running unit tests..."
    COMMAND $(MAKE) run_unit_tests
    COMMAND echo "3. Running integration tests..."
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/integration/test_enhanced.sh
    COMMAND echo "4. Running benchmarks..."
    COMMAND $(MAKE) run_benchmarks
    DEPENDS unit_tests benchmark_coverage
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
)