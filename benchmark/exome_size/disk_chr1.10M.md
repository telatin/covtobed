| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.10M ../TE/chr1.10M.bam` | 9.662 ± 0.369 | 9.314 | 10.307 | 1.44 ± 0.08 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.10M ../TE/chr1.10M.bam` | 8.801 ± 0.341 | 8.518 | 9.381 | 1.31 ± 0.07 |
| `covtobed ../TE/chr1.10M.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.10M.bed` | 6.713 ± 0.254 | 6.271 | 7.011 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.10M.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.10M.bed` | 42.749 ± 0.590 | 41.626 | 43.281 | 6.37 ± 0.26 |
