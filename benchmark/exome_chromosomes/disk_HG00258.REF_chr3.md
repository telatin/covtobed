| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr3 ../TE/HG00258.REF_chr3.bam` | 25.565 ± 0.441 | 24.842 | 26.169 | 1.32 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr3 ../TE/HG00258.REF_chr3.bam` | 19.360 ± 0.567 | 18.668 | 20.185 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr3.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr3.bed` | 89.252 ± 2.104 | 85.553 | 91.172 | 4.61 ± 0.17 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr3.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr3.bed` | 161.720 ± 3.985 | 157.259 | 165.958 | 8.35 ± 0.32 |
