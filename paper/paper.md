---
title: 'covtobed: a simple and fast tool to extract coverage tracks from BAM files'
tags:
  - bedtools
  - bamtools
  - genomics
  - bioinformatics
  - target enrichment
  - sequence coverage
authors:
  - name: Giovanni Birolo
    affiliation: 1
  - name: Andrea Telatin
    affiliation: 2
    orcid: 0000-0001-7619-281X
affiliations:
  - name: Dept. Medical Sciences, University of Turin, ITALY
    index: 1
  - name: Gut Microbes and Health Programme, Quadram Institute Bioscience, Norwich, UK
    index: 2
date: "10 December 2019"
bibliography: paper.bib

---


# Summary

A common task in bioinformatics is the mapping of DNA sequencing reads (produced by "next generation sequencing" experiments) against a reference genome. 
The output of the alignment is commonly encoded in a BAM file [@samformat]. 
For several applications of DNA sequencing it is useful to extract the *depth of coverage* [@coverage] 
at specific positions in the BAM file, 
encoding the output in the standard BED format [@bedtools].

Here we describe *covtobed*, a C++ program designed to extract the depth of coverage per position from a sorted BAM file, 
optionally specifying a range of coverage of interest and a minimum length for the features to be printed in the output BED file. 
Parsing of BAM files is performed using `libbamtools` [@bamtools]. 

The design has been inspired by the UNIX programming philosophy [@phylosophy], and thus *covtobed* performs a single task and supports input and output streams.

# Availability and Installation

*covtobed* is distributed with MIT licence and available from the [GitHub repository](https://github.com/telatin/covtobed), 
and can be easily installed via Miniconda from the "bioconda" channel 
(*i. e.* `conda install -c bioconda covtobed`).

The tool is also available as a Docker image downloadable from [Docker Hub](https://hub.docker.com/r/andreatelatin/covtobed) 
(*i. e.* `docker pull andreatelatin/covtobed`) or as a Singularity image.

# Code (structure and dependencies)

The code is object oriented, including an *Input* class handling reading, parsing and filtering of alignments and an *Output* class handling coverage filtering and writing in different formats.
The main algorithm is based on a *priority_queue* from the standard library and is both fast and memory efficient.

*covtobed* relies on  `libbamtools` [@bamtools] for BAM file parsing, and `cpp-optparse` [@opt] for command line option parsing.

# Example application

When performing target enrichment experiments (where the aim is to sequence a set of selected regions of a genome), it's important to detect a lack of coverage or unsufficient coverage (*i.e.* is lower than `THRESHOLD` in the target regions). This can be calculated intersecting (with *bedtools*, [@bedtools]) the target BED file with the output of *covtobed*. 

The tool has been used, for example, in the setup of a *target enrichment* panel targeting 71 human genes [@poloni], in order to detect uncovered regions.

While a tool exists – called *mosdepth* [@mosdepth] – to perform a coverage analysis, *covtobed* was designed with the ability to quickly extract regions between used defined coverage intervals and, more importantly, with streaming from standard input and to standard output, that Mosdepth doesn't support.

# Performance

*covtobed* is faster than its alternatives. In particular, *bedtools* is approximately 20-30 times slower, and *mosdepth* is approximately 30-60 times slower (but will produce a compressed and indexed output). We report the measured time using a BAM file from a target enrichment sequencing panel [@dataset]:

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| **mosdepth** | 66.626 ± 0.747 | 65.366 | 67.531 | 64.62 ± 1.52 |
| **covtobed** | 1.031 ± 0.021 | 1.006 | 1.070 | 1.00 |
| **bedtools** genomecov | 34.076 ± 1.644 | 33.146 | 37.360 | 33.05 ± 1.73 |

The scripts to perform the benchmark are available in the [github repository](https://github.com/telatin/covtobed/tree/master/benchmark).

# Acknowledgements

The authors gratefully acknowledge the support of the Biotechnology and Biological Sciences Research Council (BBSRC); 
this research was partly supported by the BBSRC Institute Strategic Programme Gut Microbes and Health BB/R012490/1 and its constituent project BBS/E/F/000PR10353. Analyses and benchmark performed using the MRC CLIMB cloud computing environment supported by grant MR/L015080/1.

# References

