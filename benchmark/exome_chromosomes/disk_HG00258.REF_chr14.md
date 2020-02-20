| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr14 ../TE/HG00258.REF_chr14.bam` | 13.538 ± 0.337 | 13.289 | 14.196 | 1.31 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr14 ../TE/HG00258.REF_chr14.bam` | 10.299 ± 0.242 | 9.906 | 10.575 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr14.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr14.bed` | 43.850 ± 1.902 | 41.628 | 47.358 | 4.26 ± 0.21 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr14.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr14.bed` | 102.685 ± 0.713 | 101.721 | 103.382 | 9.97 ± 0.24 |
