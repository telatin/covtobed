| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| ` mosdepth "/local/giovanni/covtobed/benchmark"/mosd2_example2 /local/giovanni/covtobed/benchmark/../TE/example2.bam` | 68.555 ± 0.765 | 67.672 | 68.999 | 33.57 ± 0.49 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/covt2_example2.bed` | 2.042 ± 0.019 | 2.024 | 2.061 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/bedt2_example2.bed` | 35.986 ± 0.167 | 35.844 | 36.169 | 17.62 ± 0.18 |
