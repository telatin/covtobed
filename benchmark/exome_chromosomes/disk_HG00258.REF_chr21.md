| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr21 ../TE/HG00258.REF_chr21.bam` | 6.511 ± 0.156 | 6.384 | 6.809 | 1.29 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr21 ../TE/HG00258.REF_chr21.bam` | 5.064 ± 0.108 | 4.965 | 5.248 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr21.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr21.bed` | 18.764 ± 0.633 | 18.115 | 19.941 | 3.71 ± 0.15 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr21.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr21.bed` | 62.386 ± 1.472 | 60.422 | 63.909 | 12.32 ± 0.39 |
