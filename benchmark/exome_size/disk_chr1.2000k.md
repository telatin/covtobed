| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.2000k ../TE/chr1.2000k.bam` | 8.392 ± 0.191 | 8.205 | 8.608 | 3.72 ± 0.13 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.2000k ../TE/chr1.2000k.bam` | 8.542 ± 0.404 | 7.834 | 9.009 | 3.79 ± 0.21 |
| `covtobed ../TE/chr1.2000k.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.2000k.bed` | 2.253 ± 0.062 | 2.176 | 2.349 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.2000k.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.2000k.bed` | 37.702 ± 0.748 | 36.839 | 38.466 | 16.73 ± 0.57 |
