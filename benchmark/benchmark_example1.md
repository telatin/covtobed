| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 974.6 ± 28.9 | 944.9 | 1020.3 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 31595.8 ± 808.6 | 30671.1 | 32461.1 | 32.42 ± 1.27 |
