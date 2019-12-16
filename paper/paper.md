---
title: 'covtobed: a fast and simple tool to extract coverage tracks from BAM files'
tags:
  - bedtools
  - bamtools
  - genomics
  - bioinformatics
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

A common task in bioinformatics is the alignment of DNA sequencing reads (produced by *next generation sequencing* experiments) against a reference genome. The output of the alignment is commonly encoded in a BAM file [@samformat], and  for several applications of DNA sequencing it is useful to extract the **depth of coverage**[@coverage] at specific positions in the BAM file, encoding the output in the standard BED format[@bedtools].

Here we describe `covtobed`, a C++ program designed to extract the depth of coverage per position from a sorted BAM file, eventually specifying a range of coverage of interest and a minimum length for the features printed in the output BED file. The fast parsing of BAM files is performed using `libbamtools` [@bamtools]. The implementation has been inspired by the UNIX programming phylosophy, and thus `covtobed` performs a limited task and supports input and output streams.

# Availability and Installation

`covtobed` is distributed with MIT licence and available from the [GitHub repository](https://github.com/telatin/covtobed), and can be easily installed via Miniconda from the "bioconda" channel (*i.e.* `conda install -c bioconda covtobed`).

The tool is also available as a Docker image downloadable from [Docker Hub](https://hub.docker.com/r/andreatelatin/covtobed) or as a Singularity image from [Zenodo](https://zenodo.org/record/1063493).

# Example application

When performing target enrichment experiments (where the aim is to sequence a set of selected regions of a genome), it's important to detect a lack of coverage or unsufficient coverage (*i.e.* is lower than `THRESHOLD` in the target regions. This can be calculated intersecting (with `bedtools`[@bedtools]) the target BED file with the output of `covtobed --max-cov THRESHOLD alignment.bam`.

While a tool exists to perform a comprehensive analysis of the coverage (`mosdepth`[@mosdepth]),  `covtobed` was designed with the ability to quickly extract regions between used defined coverage intervals.

# Acknowledgements

The initial development of the tool has been supported by the Seventh Framework Program (600932 MD-Paedigree).

The authors gratefully acknowledge the support of the Biotechnology and Biological Sciences Research Council (BBSRC); this research was partly supported by the BBSRC Institute Strategic Programme Gut Microbes and Health BB/R012490/1 and its constituent project BBS/E/F/000PR10353.


# References
