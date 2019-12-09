# covtobed

Read one (or more) alignment files (sorted BAM) and prints a BED with the coverage.

Synopsis:
```
Usage: scov [options] [BAM]...

Computes coverage from alignments

Options:
  -h, --help            show this help message and exit
  --physical-coverage   compute physical coverage (needs paired alignments in input)
  -q MIN, --min-mapq=MIN
                        skip alignments whose mapping quality is less than MIN
                        (default: 0)
  -m MINCOV, --min-cov=MINCOV
                        print BED feature only if the coverage is bigger than
                        (or equal to) MINCOV (default: 0)
  -x MAXCOV, --max-cov=MAXCOV
                        print BED feature only if the coverage is lower than
                        MAXCOV (default: 100000)
  --output-strands      outputs coverage and stats separately for each strand
  --format=CHOICE       output format
  
```

## Install

To install with Miniconda:

```
conda install -c bioconda covtobed
```
## Requirements and compiling

This tool requires **libbamtools** and **zlib**.

To manually compile:
```
c++ -std=c++11 *.cpp -I/path/to/bamtools/ -L${HOME}/path/to/lib/ -lbamtools -o covtobed
```

