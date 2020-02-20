| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr16 ../TE/HG00258.REF_chr16.bam` | 14.139 ± 0.135 | 13.943 | 14.278 | 1.39 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr16 ../TE/HG00258.REF_chr16.bam` | 10.181 ± 0.249 | 9.938 | 10.649 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr16.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr16.bed` | 47.331 ± 1.437 | 46.201 | 50.049 | 4.65 ± 0.18 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr16.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr16.bed` | 108.883 ± 1.818 | 106.251 | 110.937 | 10.69 ± 0.32 |
