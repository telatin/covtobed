| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.5000k ../TE/chr1.5000k.bam` | 9.223 ± 0.318 | 8.922 | 9.787 | 2.34 ± 0.11 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.5000k ../TE/chr1.5000k.bam` | 8.998 ± 0.374 | 8.310 | 9.352 | 2.29 ± 0.12 |
| `covtobed ../TE/chr1.5000k.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.5000k.bed` | 3.936 ± 0.123 | 3.812 | 4.168 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.5000k.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.5000k.bed` | 40.330 ± 0.773 | 38.870 | 41.037 | 10.25 ± 0.38 |
