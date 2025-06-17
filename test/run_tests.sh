#!/bin/bash

# Convenience script for running the enhanced test suite

set -eou pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Parse command line arguments
MODE="all"
BUILD_DIR="build"

while [[ $# -gt 0 ]]; do
    case $1 in
        -u|--unit)
            MODE="unit"
            shift
            ;;
        -i|--integration)
            MODE="integration"
            shift
            ;;
        -b|--benchmark)
            MODE="benchmark"
            shift
            ;;
        -a|--all)
            MODE="all"
            shift
            ;;
        --build-dir)
            BUILD_DIR="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  -u, --unit         Run unit tests only"
            echo "  -i, --integration  Run integration tests only"
            echo "  -b, --benchmark    Run benchmarks only"
            echo "  -a, --all          Run all tests (default)"
            echo "  --build-dir DIR    Specify build directory (default: build)"
            echo "  -h, --help         Show this help"
            exit 0
            ;;
        *)
            print_error "Unknown option: $1"
            exit 1
            ;;
    esac
done

print_header "Covtobed Enhanced Test Suite"

# Check if covtobed binary exists
if [ ! -f "../covtobed" ]; then
    print_warning "covtobed binary not found. Building..."
    cd ..
    c++ -std=c++11 *.cpp -I/usr/include/bamtools -lbamtools -o covtobed -lz
    cd test
fi

# Create build directory
if [ ! -d "$BUILD_DIR" ]; then
    print_status "Creating build directory: $BUILD_DIR"
    mkdir -p "$BUILD_DIR"
fi

cd "$BUILD_DIR"

# Configure with CMake
if [ ! -f "Makefile" ]; then
    print_status "Configuring with CMake..."
    cmake ..
fi

# Build test suite
print_status "Building test suite..."
make -j$(nproc)

# Run tests based on mode
case $MODE in
    "unit")
        print_header "Running Unit Tests"
        make run_unit_tests
        ;;
    "integration")
        print_header "Running Integration Tests"
        cd ..
        ./integration/test_enhanced.sh
        ;;
    "benchmark")
        print_header "Running Benchmarks"
        make run_benchmarks
        ;;
    "all")
        print_header "Running Comprehensive Test Suite"
        
        # Generate test data
        print_status "Generating synthetic test data..."
        make generate_test_data
        
        # Run unit tests
        print_status "Running unit tests..."
        make run_unit_tests
        
        # Run integration tests
        print_status "Running integration tests..."
        cd ..
        ./integration/test_enhanced.sh
        cd "$BUILD_DIR"
        
        # Run benchmarks
        print_status "Running benchmarks..."
        make run_benchmarks
        
        print_header "All Tests Completed Successfully!"
        ;;
esac

cd "$SCRIPT_DIR"
print_status "Test suite execution completed."