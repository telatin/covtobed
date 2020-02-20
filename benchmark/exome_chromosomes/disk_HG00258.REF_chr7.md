| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr7 ../TE/HG00258.REF_chr7.bam` | 20.934 ± 0.442 | 20.504 | 21.773 | 1.34 ± 0.03 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr7 ../TE/HG00258.REF_chr7.bam` | 15.579 ± 0.180 | 15.310 | 15.742 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr7.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr7.bed` | 74.747 ± 1.936 | 71.174 | 76.559 | 4.80 ± 0.14 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr7.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr7.bed` | 139.744 ± 2.427 | 136.704 | 143.761 | 8.97 ± 0.19 |
