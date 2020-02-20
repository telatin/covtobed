| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr18 ../TE/HG00258.REF_chr18.bam` | 10.144 ± 0.177 | 9.860 | 10.310 | 1.24 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr18 ../TE/HG00258.REF_chr18.bam` | 8.154 ± 0.218 | 7.885 | 8.410 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr18.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr18.bed` | 32.654 ± 1.230 | 31.418 | 34.061 | 4.00 ± 0.18 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr18.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr18.bed` | 84.882 ± 1.210 | 83.370 | 86.689 | 10.41 ± 0.32 |
