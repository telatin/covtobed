# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`covtobed` is a C++ bioinformatics tool that generates BED coverage tracks from sorted BAM alignment files. It processes sequence alignment data to compute coverage depth and output regions in BED format, with support for physical coverage, strand-specific analysis, and various filtering options.

## Build System and Dependencies

### Manual Compilation
The project uses direct C++ compilation without makefiles:

```bash
# Basic compilation
c++ -std=c++11 *.cpp -I/path/to/bamtools/ -L${HOME}/path/to/lib/ -lbamtools -o covtobed

# Ubuntu with system packages
g++ -std=c++11 *.cpp -I/usr/include/bamtools /usr/lib/x86_64-linux-gnu/libbamtools.a -o covtobed -lz

# macOS with conda
c++ -std=c++11 *.cpp -I"${HOME}"/miniconda3/include/bamtools/ -L"${HOME}"/miniconda3/lib/ "${HOME}"/miniconda3/lib/libbamtools.a -o covtobed -lz
```

### Dependencies
- **libbamtools**: Required for BAM file handling
- **zlib**: For compression support
- **C++11**: Minimum standard required

### Platform-Specific Build Scripts
- `binaries/build_ubuntu.sh`: Automated Ubuntu build with static linking
- `binaries/build_osx.sh`: Automated macOS build using conda dependencies

## Testing

### Test Suite
Run the comprehensive test suite:
```bash
bash test/test.sh
```

The test script:
- Uses the `covtobed` binary in the project root
- Falls back to pre-compiled binaries if main binary doesn't exist
- Tests various functionality: sorted/unsorted BAM handling, coverage filtering, physical coverage, strand analysis, output formats
- Compares output against expected results in `test/output.bed` and `test/mock.bed`
- Sets `COVTOBED_QUIET=1` environment variable to suppress startup messages

### Test Data
- `test/demo.bam`: Main test file for standard coverage analysis
- `test/mp.bam`: Mate-pair BAM for physical coverage testing
- `test/mock.bam`: Synthetic BAM with known coverage values
- `test/stranded.bam`: Test file for strand-specific analysis
- Additional test files for edge cases (duplicates, filtering, etc.)

## Code Architecture

### Core Components
- **base.cpp**: Main program logic, coverage calculation engine
  - Uses BamTools library for BAM I/O
  - Implements coverage tracking with priority queues
  - Handles strand-specific and physical coverage modes
- **OptionParser.cpp/h**: Command-line argument parsing
- **interval.h**: Data structures for genomic intervals
- **covtobed**: Main executable (compiled binary)

### Key Classes and Structures
- `Input`: Handles BAM file reading and alignment filtering
- `CovEnd`: Manages alignment end positions in priority queue
- `DepthType`: Type definition for coverage depth (uint32_t)
- Coverage calculation uses STL priority queues for efficient processing

### Debug Mode
Debug output can be enabled by changing `#define debug if(false)` to `#define debug if(true)` in base.cpp before compilation.

## CI/CD

The project uses GitHub Actions (`.github/workflows/c-cpp.yml`):
- Builds on Ubuntu 18.04
- Installs libbamtools-dev system package
- Compiles with both g++ and clang++
- Runs basic functionality tests and full test suite

## Development Workflow

1. **Building**: Use platform-appropriate build script or manual compilation command
2. **Testing**: Always run `bash test/test.sh` before submitting changes
3. **Debugging**: Enable debug mode in base.cpp for detailed output
4. **Pull Requests**: Follow GitHub Flow - fork, branch, PR against master

## Environment Variables

- `COVTOBED_QUIET=1`: Suppresses startup message when reading from STDIN

## Key Features to Understand

- **Sorted BAM Requirement**: Tool requires sorted BAM input, will error on unsorted files
- **Streaming Support**: Can read from STDIN for pipeline integration
- **Physical Coverage**: Special mode for paired-end/mate-pair libraries
- **Strand Analysis**: Can output strand-specific coverage information
- **Multiple Output Formats**: BED (default) and counts formats supported
- **Filtering Options**: By mapping quality, coverage thresholds, alignment validity