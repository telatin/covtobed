| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr8 ../TE/HG00258.REF_chr8.bam` | 17.817 ± 0.273 | 17.434 | 18.138 | 1.29 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr8 ../TE/HG00258.REF_chr8.bam` | 13.808 ± 0.391 | 13.196 | 14.238 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr8.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr8.bed` | 62.683 ± 2.771 | 57.827 | 65.977 | 4.54 ± 0.24 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr8.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr8.bed` | 124.332 ± 2.761 | 122.286 | 129.520 | 9.00 ± 0.32 |
