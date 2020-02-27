| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.150M ../TE/chr1.150M.bam` | 24.784 ± 0.668 | 23.808 | 25.531 | 1.35 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.150M ../TE/chr1.150M.bam` | 18.390 ± 0.483 | 17.728 | 18.970 | 1.00 |
| `covtobed ../TE/chr1.150M.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.150M.bed` | 76.797 ± 1.990 | 73.880 | 79.782 | 4.18 ± 0.15 |
| `bedtools genomecov -bga -ibam ../TE/chr1.150M.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.150M.bed` | 144.268 ± 4.728 | 139.624 | 151.378 | 7.85 ± 0.33 |
