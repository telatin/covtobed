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

The sequencing of the human genome an

Source Code: https://github.com/telatin/covtobed.
Contact: andrea.telatin@quadram.ac.uk


# Summary

A common task in bioinformatics is the alignment of DNA sequencing reads against a reference genome, commonly encoded in a BAM file [@Li-sam]. For several applications of DNA sequencing it is useful to have a report on the coverage [].
The sequencing of the human genome and subsequent advances in DNA sequencing technology have transformed modern biological research by producing data sets of ever-increasing size and complexity. While these technologies have led to breakthroughs in genetics research, the incredible throughput and breadth of the resulting data have spurred a reliance on computational tools and programming languages to interpret the results. In 2010, Quinlan et al. developed bedtools, a powerful suite of command-line tools for ‘genome arithmetic’ that has become one of the most widely used and indispensable tools for genomic data analysis [@Quinlan2010-ak]. A year later, `pybedtools` extended the features of bedtools to the python programming language [@Dale2011-ja]. During that same time period, the use of the programming language R—with a rich trove of libraries for statistical analysis and data visualization—skyrocketed in the biological sciences [@Tippmann2015-yz]. While several R packages have been developed for bedtools-like genome analysis, their usage and functionality differ significantly from that of bedtools [@Riemondy2017-vz; @Lawrence2013-px]. These differences make them more difficult to use for those who are already familiar with bedtools behavior and lacks some of the capabilities of bedtools.

Here we describe `covtobed`, an R package that allows seamless integration of bedtools functions into the R programming environment. `covtobed` functionality, inputs, outputs, and documentation perfectly replicate those found in the command-line version of bedtools and offer new features for improved usability within the R environment.



# Discussion

`covtobed` provides seamless integration of the bedtools software suite into the R programming environment. The package was designed to be as user-friendly as possible and should be intuitive for those already familiar with the bedtools command-line tool. The ability to handle multiple data types, the forward and backward compatibility, and the included unit tests ensure stability and ease of use. The harmonious combination of these two powerful analytical platforms should make `covtobed` a valuable and widely used tool for genomic analysis.



# References

