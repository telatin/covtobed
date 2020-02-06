| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `samtools index /local/giovanni/covtobed/benchmark/../TE/example1.bam && mosdepth "/local/giovanni/covtobed/benchmark"/mosd_example1 /local/giovanni/covtobed/benchmark/../TE/example1.bam && gunzip "/local/giovanni/covtobed/benchmark"/mosd_example1.per-base.bed.gz` | 70.926 ± 2.220 | 69.095 | 73.395 | 65.69 ± 2.14 |
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example1.bam > "/local/giovanni/covtobed/benchmark"/covt_example1.bed` | 1.080 ± 0.010 | 1.071 | 1.090 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example1.bam > "/local/giovanni/covtobed/benchmark"/bedt_example1.bed` | 36.506 ± 2.616 | 34.523 | 39.471 | 33.81 ± 2.44 |
