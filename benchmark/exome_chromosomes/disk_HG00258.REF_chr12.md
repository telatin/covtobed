| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr12 ../TE/HG00258.REF_chr12.bam` | 19.543 ± 0.410 | 19.024 | 20.173 | 1.39 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr12 ../TE/HG00258.REF_chr12.bam` | 14.081 ± 0.288 | 13.858 | 14.509 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr12.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr12.bed` | 67.944 ± 1.702 | 65.743 | 70.557 | 4.83 ± 0.16 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr12.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr12.bed` | 138.625 ± 2.803 | 135.714 | 143.055 | 9.84 ± 0.28 |
