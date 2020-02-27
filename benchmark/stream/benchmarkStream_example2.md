| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `samtools depth -a ../TE/example2.bam > /dev/null` | 606.278 ± 1.351 | 604.367 | 607.973 | 308.35 ± 6.71 |
| `samtools depth ../TE/example2.bam > /dev/null` | 8.779 ± 0.152 | 8.609 | 9.002 | 4.47 ± 0.12 |
| `covtobed ../TE/example2.bam > /dev/null` | 1.966 ± 0.043 | 1.917 | 2.011 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/example2.bam > /dev/null` | 37.296 ± 0.202 | 37.094 | 37.571 | 18.97 ± 0.42 |
