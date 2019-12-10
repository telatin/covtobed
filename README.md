# covtobed - a tool to generate BED coverage tracks from BAM files

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/covtobed/README.html)
![License](https://anaconda.org/bioconda/covtobed/badges/license.svg)
![Platforms](https://anaconda.org/bioconda/covtobed/badges/platforms.svg)
![Downloads](https://anaconda.org/bioconda/covtobed/badges/downloads.svg)

Read one (or more) alignment files (sorted BAM) and prints a BED with the coverage. It will join consecutive bases with the same coverage, and can be used to only print a BED file with the regions having a specific coverage range.


![covtobed example](img/coverage_bam_to_bed.png)



Synopsis:
```
Usage: covtobed [options] [BAM]...

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
## Example

Command:
```
covtobed -m 0 -x 5 test/demo.bam
```
Output:
```
NC_001416.1	1112	1120	4
NC_001416.1	1120	1999	0
NC_001416.1	2653	2674	4
NC_001416.1	2674	2686	3
NC_001416.1	2686	2762	2
NC_001416.1	2762	2771	1
NC_001416.1	2771	2800	0
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
## Acknowledgements

This tools uses [libbamtools](https://github.com/pezmaster31/bamtools) by Derek Barnett, Erik Garrison, Gabor Marth and Michael Stromberg, and [cpp-optparse](https://github.com/weisslj/cpp-optparse) by Johannes Weißl. Both tools and this program are released with MIT license.

## Authors

Giovanni Birolo ([@gbirolo](https://github.com/gbirolo)) and Andrea Telatin ([@telatin](https://github.com/telatin))
