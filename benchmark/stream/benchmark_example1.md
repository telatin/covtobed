| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `samtools depth -a /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 600.603 ± 3.205 | 596.251 | 604.440 | 573.92 ± 26.79 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 1.046 ± 0.049 | 0.982 | 1.096 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 34.717 ± 0.952 | 33.325 | 36.090 | 33.17 ± 1.79 |
