cmake_minimum_required(VERSION 3.10)
project(covtobed_unit_tests)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find required packages
find_package(PkgConfig REQUIRED)
pkg_check_modules(BAMTOOLS REQUIRED bamtools-1)
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    find_package(Catch2 2 REQUIRED)
endif()

# Include directories
include_directories(${BAMTOOLS_INCLUDE_DIRS})
include_directories(../..)

# Create unit test executable
add_executable(unit_tests
    test_main.cpp
    test_coverage.cpp
    test_intervals.cpp
    test_output.cpp
)

# Link libraries
if(TARGET Catch2::Catch2WithMain)
    target_link_libraries(unit_tests Catch2::Catch2WithMain)
else()
    target_link_libraries(unit_tests Catch2::Catch2)
endif()

target_link_libraries(unit_tests 
    ${BAMTOOLS_LIBRARIES}
    z
)

# Add compiler flags
target_compile_options(unit_tests PRIVATE ${BAMTOOLS_CFLAGS_OTHER})

# Enable testing
enable_testing()
add_test(NAME unit_tests COMMAND unit_tests)

# Custom target to run tests
add_custom_target(run_unit_tests
    COMMAND unit_tests
    DEPENDS unit_tests
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)