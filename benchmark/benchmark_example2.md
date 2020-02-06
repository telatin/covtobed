| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `samtools index /local/giovanni/covtobed/benchmark/../TE/example2.bam && mosdepth "/local/giovanni/covtobed/benchmark"/mosd_example2 /local/giovanni/covtobed/benchmark/../TE/example2.bam && gunzip "/local/giovanni/covtobed/benchmark"/mosd_example2.per-base.bed.gz` | 70.060 ± 1.545 | 69.001 | 71.833 | 34.67 ± 0.79 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/covt_example2.bed` | 2.021 ± 0.012 | 2.013 | 2.034 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example2.bam > "/local/giovanni/covtobed/benchmark"/bedt_example2.bed` | 37.959 ± 1.524 | 36.790 | 39.682 | 18.79 ± 0.76 |
