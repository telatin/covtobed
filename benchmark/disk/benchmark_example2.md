| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_example2 /local/giovanni/covtobed/benchmark/../TE/example2.bam` | 67.378 ± 0.952 | 65.847 | 68.399 | 35.58 ± 0.64 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_example2 /local/giovanni/covtobed/benchmark/../TE/example2.bam` | 68.747 ± 1.801 | 67.027 | 71.854 | 36.30 ± 1.03 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/covt2_example2.bed` | 1.894 ± 0.021 | 1.871 | 1.919 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/bedt2_example2.bed` | 37.326 ± 1.609 | 35.432 | 40.154 | 19.71 ± 0.88 |
