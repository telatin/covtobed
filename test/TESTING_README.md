# Enhanced Testing Suite for covtobed

This directory contains a comprehensive testing framework for the covtobed bioinformatics tool, 
including unit tests, integration tests, benchmarks, and test data generators.

## Quick Start

### Ubuntu dependencies

```bash
# Core dependencies
sudo apt install -y build-essential cmake pkg-config libbamtools-dev libjsoncpp-dev zlib1g-dev

# Testing frameworks
sudo apt install -y catch2 libbenchmark-dev samtools valgrind bc
```
### Run All Tests
```bash
cd test
bash run_tests.sh --all
```

### Run Specific Test Types
```bash
# Unit tests only
bash run_tests.sh --unit

# Integration tests only  
bash run_tests.sh --integration

# Benchmarks only
bash run_tests.sh --benchmark
```

### Using CMake

```bash
cd test
mkdir build && cd build
cmake ..
make

# Run specific targets
make run_unit_tests
make run_benchmarks
make generate_test_data
```

## 📁 Directory Structure

```
test/
├── unit/                     # Unit tests with Catch2
│   ├── test_coverage.cpp     # Tests for Coverage class
│   ├── test_intervals.cpp    # Tests for interval operations
│   ├── test_output.cpp       # Tests for output formatting
│   └── CMakeLists.txt
├── integration/              # Enhanced integration tests
│   ├── test_enhanced.sh      # Comprehensive integration tests
│   └── test_helpers.sh       # Test helper functions
├── benchmark/                # Performance benchmarks
│   ├── benchmark_coverage.cpp# Core algorithm benchmarks
│   └── CMakeLists.txt
├── utils/                    # Test utilities
│   ├── simple_test_generator.sh # SAM/BAM test data generator
│   └── synthetic_bam_generator.cpp # (C++ generator - disabled)
├── data/                     # Test data
│   ├── synthetic/            # Generated test files
│   └── reference_outputs/    # Expected outputs
├── run_tests.sh              # Main test runner script
└── CMakeLists.txt            # Master build configuration
```

## 🧪 Unit Tests

### Coverage
- **Framework**: Catch2 3.x
- **Files**: `test_coverage.cpp`, `test_intervals.cpp`, `test_output.cpp`
- **Focus**: Core classes (Coverage, CovEnd, Interval)

### Running Unit Tests
```bash
cd test/build
make run_unit_tests
# or
./unit/unit_tests
```
 
## 🔗 Integration Tests

### Enhanced Integration Testing
- **Framework**: Bash with colored output and detailed reporting
- **File**: `test_enhanced.sh`
- **Features**: 
  - Comprehensive error reporting
  - Performance measurement
  - Memory usage testing (with valgrind)
  - Edge case validation
 ### Helper Functions
The `test_helpers.sh` provides:
- `run_test()`: Execute tests with error capture
- `assert_*()`: Various assertion functions
- `measure_time()`: Performance measurement
- `print_status()`: Colored output formatting

## 📊 Benchmarks

### Performance Testing
- **Framework**: Google Benchmark
- **File**: `benchmark_coverage.cpp`
- **Focus**: Core algorithm performance

### Benchmark Categories
1. **Coverage operations**: increment/decrement performance
2. **Priority queue**: operations with varying data sizes
3. **Interval operations**: intersection, difference calculations
4. **Memory allocation**: allocation patterns and performance
5. **Coverage simulation**: end-to-end algorithm performance

### Running Benchmarks
```bash
cd test/build
make run_benchmarks
# or
./benchmark/benchmark_coverage
```

### Sample Output
```
BM_CoverageIncrement                        1.38 ns         1.38 ns    506618768
BM_PriorityQueueOperations/100              4330 ns         4328 ns       159753
BM_CoverageCalculationSimulation/1000        2.5 ms          2.5 ms          280
```

## 🗃️ Test Data Generation

### Synthetic Data Generator
- **Script**: `simple_test_generator.sh`
- **Formats**: SAM (with optional BAM conversion via samtools)

### Generation Types
```bash
# Specific coverage pattern
./simple_test_generator.sh output.sam coverage 100 200 10

# Stranded coverage data
./simple_test_generator.sh output.sam stranded

# Single read file
./simple_test_generator.sh output.sam single
```

### Automated Data Generation
```bash
make generate_test_data
```

## 🛠️ Dependencies

### Required Tools
- **catch2**: Unit testing framework
- **libbenchmark-dev**: Performance benchmarking
- **pkg-config**: Build configuration
- **libjsoncpp-dev**: JSON support for bamtools
- **libbamtools-dev**: BAM file handling

### Optional Tools
- **samtools**: SAM/BAM conversion for test data
- **valgrind**: Memory analysis (for integration tests)
- **bc**: Floating point math in shell scripts

### Installation
```bash
sudo apt update
sudo apt install -y catch2 libbenchmark-dev cmake pkg-config libjsoncpp-dev libbamtools-dev
```

## 📈 Performance Monitoring

### Benchmark Results Storage
```bash
# Save benchmark results to JSON
./benchmark/benchmark_coverage --benchmark_format=json --benchmark_out=results.json

# Compare with previous results
# (requires benchmark comparison tools)
```

### Integration Test Performance
The enhanced integration tests include basic performance monitoring:
- Execution time measurement
- Memory usage tracking (with valgrind)
- Performance regression detection

## 🐛 Debugging Tests

### Unit Test Debugging
```bash
# Run specific test cases
./unit/unit_tests "Coverage class"

# Run with verbose output
./unit/unit_tests -v

# List all test cases
./unit/unit_tests --list-tests
```

### Integration Test Debugging
```bash
# Run single test
TESTS_RUN=0 TESTS_PASSED=0 TESTS_FAILED=0
source test/integration/test_helpers.sh
run_test "test_name" "command"
```

## 🔄 Continuous Testing

### Local Development Workflow
1. Make code changes
2. Run unit tests: `make run_unit_tests`
3. Run integration tests: `bash test/integration/test_enhanced.sh`
4. Run benchmarks: `make run_benchmarks`
5. Check performance regressions

### Test Coverage Analysis
```bash
# Compile with coverage flags (add to CMakeLists.txt)
# -fprofile-arcs -ftest-coverage

# Generate coverage report
# gcov test/unit/*.cpp
# lcov --capture --directory . --output-file coverage.info
```

## ✅ Test Quality Metrics

### Current Status
- **Unit Tests**: ✅ Framework implemented, basic coverage classes tested
- **Integration Tests**: ✅ 18 comprehensive test cases
- **Benchmarks**: ✅ 8 performance benchmarks covering core algorithms
- **Test Data**: ✅ Synthetic data generation for various scenarios

### Future Improvements
- [ ] Increase unit test coverage to 80%+
- [ ] Add property-based testing for interval operations
- [ ] Implement performance regression detection
- [ ] Add fuzzing tests for robustness
- [ ] Cross-platform testing automation

## 📚 Contributing

### Adding New Tests

1. **Unit Tests**: Add to appropriate `test_*.cpp` file in `test/unit/`
2. **Integration Tests**: Add to `test_enhanced.sh` using helper functions
3. **Benchmarks**: Add to `benchmark_coverage.cpp` with appropriate naming

### Test Naming Conventions
- Unit tests: `TEST_CASE("Description", "[tag]")`
- Integration tests: `run_test "test_name" "command"`
- Benchmarks: `BM_OperationName` or `BM_OperationName/Parameter`

### Code Style
- Follow existing patterns in test files
- Use descriptive test names and clear assertions
- Include both positive and negative test cases
- Add comments for complex test scenarios

This enhanced testing suite provides a solid foundation for maintaining and improving the covtobed codebase quality.