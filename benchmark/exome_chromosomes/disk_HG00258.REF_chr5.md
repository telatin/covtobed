| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr5 ../TE/HG00258.REF_chr5.bam` | 22.886 ± 0.410 | 22.464 | 23.470 | 1.34 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr5 ../TE/HG00258.REF_chr5.bam` | 17.050 ± 0.357 | 16.555 | 17.625 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr5.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr5.bed` | 78.351 ± 2.589 | 74.010 | 80.779 | 4.60 ± 0.18 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr5.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr5.bed` | 151.817 ± 3.474 | 146.590 | 155.722 | 8.90 ± 0.28 |
