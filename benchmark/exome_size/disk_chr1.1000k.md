| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.1000k ../TE/chr1.1000k.bam` | 8.085 ± 0.128 | 7.923 | 8.203 | 5.69 ± 0.20 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.1000k ../TE/chr1.1000k.bam` | 8.302 ± 0.289 | 8.124 | 8.861 | 5.84 ± 0.27 |
| `covtobed ../TE/chr1.1000k.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.1000k.bed` | 1.422 ± 0.044 | 1.371 | 1.499 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.1000k.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.1000k.bed` | 37.181 ± 0.534 | 36.333 | 37.857 | 26.15 ± 0.89 |
