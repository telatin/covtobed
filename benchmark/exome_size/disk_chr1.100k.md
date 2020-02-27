| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.100k ../TE/chr1.100k.bam` | 8.114 ± 0.252 | 7.816 | 8.446 | 8.41 ± 0.32 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.100k ../TE/chr1.100k.bam` | 8.294 ± 0.333 | 7.949 | 8.861 | 8.60 ± 0.40 |
| `covtobed ../TE/chr1.100k.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.100k.bed` | 0.965 ± 0.022 | 0.942 | 1.004 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.100k.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.100k.bed` | 35.306 ± 0.295 | 34.918 | 35.607 | 36.60 ± 0.89 |
