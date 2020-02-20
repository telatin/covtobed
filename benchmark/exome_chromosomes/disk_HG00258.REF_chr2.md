| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr2 ../TE/HG00258.REF_chr2.bam` | 30.473 ± 0.343 | 30.070 | 31.065 | 1.31 ± 0.02 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr2 ../TE/HG00258.REF_chr2.bam` | 23.329 ± 0.192 | 23.015 | 23.488 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr2.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr2.bed` | 111.029 ± 2.747 | 107.657 | 114.831 | 4.76 ± 0.12 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr2.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr2.bed` | 196.879 ± 5.299 | 190.250 | 204.934 | 8.44 ± 0.24 |
