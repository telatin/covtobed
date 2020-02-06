| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_example1 /local/giovanni/covtobed/benchmark/../TE/example1.bam` | 66.626 ± 0.747 | 65.366 | 67.531 | 64.62 ± 1.52 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_example1 /local/giovanni/covtobed/benchmark/../TE/example1.bam` | 66.100 ± 0.718 | 65.654 | 67.552 | 64.11 ± 1.50 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example1.bam > "/local/giovanni/covtobed/benchmark"/covt2_example1.bed` | 1.031 ± 0.021 | 1.006 | 1.070 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example1.bam > "/local/giovanni/covtobed/benchmark"/bedt2_example1.bed` | 34.076 ± 1.644 | 33.146 | 37.360 | 33.05 ± 1.73 |
