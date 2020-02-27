| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.50M ../TE/chr1.50M.bam` | 14.940 ± 0.077 | 14.827 | 15.018 | 1.20 ± 0.02 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.50M ../TE/chr1.50M.bam` | 12.464 ± 0.206 | 12.288 | 12.801 | 1.00 |
| `covtobed ../TE/chr1.50M.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.50M.bed` | 33.034 ± 0.620 | 32.566 | 34.207 | 2.65 ± 0.07 |
| `bedtools genomecov -bga -ibam ../TE/chr1.50M.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.50M.bed` | 83.922 ± 1.180 | 82.598 | 85.508 | 6.73 ± 0.15 |
