| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr15 ../TE/HG00258.REF_chr15.bam` | 13.861 ± 0.126 | 13.694 | 13.988 | 1.34 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr15 ../TE/HG00258.REF_chr15.bam` | 10.318 ± 0.325 | 9.938 | 10.824 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr15.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr15.bed` | 45.066 ± 1.360 | 43.241 | 46.574 | 4.37 ± 0.19 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr15.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr15.bed` | 103.312 ± 2.433 | 100.651 | 106.629 | 10.01 ± 0.39 |
