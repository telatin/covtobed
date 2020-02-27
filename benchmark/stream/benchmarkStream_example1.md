| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `samtools depth -a ../TE/example1.bam > /dev/null` | 600.513 ± 3.980 | 594.494 | 605.576 | 605.45 ± 6.28 |
| `samtools depth ../TE/example1.bam > /dev/null` | 2.025 ± 0.075 | 1.970 | 2.169 | 2.04 ± 0.08 |
| `covtobed ../TE/example1.bam > /dev/null` | 0.992 ± 0.008 | 0.981 | 1.001 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/example1.bam > /dev/null` | 34.376 ± 1.032 | 33.225 | 35.964 | 34.66 ± 1.08 |
