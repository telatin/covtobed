| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr20 ../TE/HG00258.REF_chr20.bam` | 10.187 ± 0.283 | 9.910 | 10.709 | 1.30 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr20 ../TE/HG00258.REF_chr20.bam` | 7.843 ± 0.180 | 7.650 | 8.126 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr20.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr20.bed` | 34.802 ± 0.810 | 33.555 | 35.792 | 4.44 ± 0.14 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr20.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr20.bed` | 83.705 ± 0.782 | 82.808 | 84.709 | 10.67 ± 0.26 |
