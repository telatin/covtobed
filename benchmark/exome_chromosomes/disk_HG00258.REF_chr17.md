| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr17 ../TE/HG00258.REF_chr17.bam` | 15.709 ± 0.425 | 15.105 | 16.127 | 1.44 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr17 ../TE/HG00258.REF_chr17.bam` | 10.940 ± 0.282 | 10.570 | 11.349 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr17.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr17.bed` | 55.456 ± 1.012 | 54.118 | 56.980 | 5.07 ± 0.16 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr17.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr17.bed` | 118.074 ± 2.926 | 114.742 | 121.235 | 10.79 ± 0.39 |
