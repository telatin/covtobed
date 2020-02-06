| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 986.2 ± 9.5 | 977.1 | 1003.5 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example1.bam > /dev/null` | 33596.2 ± 1550.4 | 31982.4 | 35700.0 | 34.07 ± 1.61 |
