| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr11 ../TE/HG00258.REF_chr11.bam` | 19.719 ± 0.419 | 19.254 | 20.390 | 1.39 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr11 ../TE/HG00258.REF_chr11.bam` | 14.165 ± 0.212 | 13.991 | 14.570 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr11.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr11.bed` | 72.303 ± 3.274 | 69.919 | 77.485 | 5.10 ± 0.24 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr11.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr11.bed` | 137.249 ± 3.279 | 131.715 | 141.087 | 9.69 ± 0.27 |
