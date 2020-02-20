| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr19 ../TE/HG00258.REF_chr19.bam` | 12.805 ± 0.166 | 12.538 | 13.048 | 1.45 ± 0.03 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr19 ../TE/HG00258.REF_chr19.bam` | 8.844 ± 0.177 | 8.553 | 9.064 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr19.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr19.bed` | 43.306 ± 0.834 | 42.432 | 44.641 | 4.90 ± 0.14 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr19.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr19.bed` | 103.669 ± 1.743 | 101.215 | 105.787 | 11.72 ± 0.31 |
