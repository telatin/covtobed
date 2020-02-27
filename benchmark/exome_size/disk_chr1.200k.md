| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.200k ../TE/chr1.200k.bam` | 8.477 ± 0.469 | 8.013 | 9.187 | 8.08 ± 0.59 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.200k ../TE/chr1.200k.bam` | 8.263 ± 0.257 | 7.966 | 8.614 | 7.88 ± 0.45 |
| `covtobed ../TE/chr1.200k.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.200k.bed` | 1.049 ± 0.051 | 0.987 | 1.114 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.200k.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.200k.bed` | 35.882 ± 1.230 | 34.409 | 37.515 | 34.21 ± 2.02 |
