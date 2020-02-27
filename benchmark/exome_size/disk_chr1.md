| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1 ../TE/chr1.bam` | 35.884 ± 0.777 | 34.759 | 36.847 | 1.40 ± 0.06 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1 ../TE/chr1.bam` | 25.692 ± 1.011 | 24.265 | 27.199 | 1.00 |
| `covtobed ../TE/chr1.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.bed` | 125.842 ± 1.983 | 123.868 | 128.292 | 4.90 ± 0.21 |
| `bedtools genomecov -bga -ibam ../TE/chr1.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.bed` | 224.625 ± 4.205 | 219.691 | 230.415 | 8.74 ± 0.38 |
