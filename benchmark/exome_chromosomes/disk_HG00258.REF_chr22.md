| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr22 ../TE/HG00258.REF_chr22.bam` | 7.509 ± 0.108 | 7.368 | 7.640 | 1.29 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr22 ../TE/HG00258.REF_chr22.bam` | 5.825 ± 0.162 | 5.635 | 6.093 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr22.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr22.bed` | 22.990 ± 0.722 | 21.728 | 23.666 | 3.95 ± 0.17 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr22.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr22.bed` | 66.403 ± 1.534 | 64.734 | 68.679 | 11.40 ± 0.41 |
