| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_example2 /local/giovanni/covtobed/benchmark/../TE/example2.bam` | 67.692 ± 1.221 | 66.288 | 68.509 | 32.89 ± 0.60 |
| `mosdepth -x -t 8 "/local/giovanni/covtobed/benchmark"/mosd2_example2 /local/giovanni/covtobed/benchmark/../TE/example2.bam` | 67.681 ± 0.954 | 66.974 | 68.766 | 32.89 ± 0.47 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/covt2_example2.bed` | 2.058 ± 0.006 | 2.054 | 2.065 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/bedt2_example2.bed` | 37.677 ± 0.962 | 36.568 | 38.281 | 18.31 ± 0.47 |
