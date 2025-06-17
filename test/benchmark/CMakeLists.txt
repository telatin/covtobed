cmake_minimum_required(VERSION 3.10)
project(covtobed_benchmarks)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)
pkg_check_modules(BAMTOOLS REQUIRED bamtools-1)
find_package(benchmark REQUIRED)

# Include directories
include_directories(${BAMTOOLS_INCLUDE_DIRS})
include_directories(../..)

# Create benchmark executable
add_executable(benchmark_coverage
    benchmark_coverage.cpp
)

# Link libraries
target_link_libraries(benchmark_coverage 
    benchmark::benchmark
    ${BAMTOOLS_LIBRARIES}
    z
    pthread
)

# Add compiler flags
target_compile_options(benchmark_coverage PRIVATE 
    ${BAMTOOLS_CFLAGS_OTHER}
    -O3
    -DNDEBUG
)

# Custom target to run benchmarks
add_custom_target(run_benchmarks
    COMMAND benchmark_coverage --benchmark_format=console
    DEPENDS benchmark_coverage
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

# Target to run benchmarks and save results
add_custom_target(benchmark_results
    COMMAND benchmark_coverage --benchmark_format=json --benchmark_out=benchmark_results.json
    DEPENDS benchmark_coverage
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)