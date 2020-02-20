| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr4 ../TE/HG00258.REF_chr4.bam` | 22.520 ± 0.421 | 21.781 | 22.984 | 1.31 ± 0.03 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr4 ../TE/HG00258.REF_chr4.bam` | 17.208 ± 0.282 | 16.806 | 17.645 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr4.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr4.bed` | 76.963 ± 3.382 | 73.592 | 82.883 | 4.47 ± 0.21 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr4.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr4.bed` | 144.057 ± 4.126 | 138.739 | 149.479 | 8.37 ± 0.28 |
