| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr9 ../TE/HG00258.REF_chr9.bam` | 17.865 ± 0.289 | 17.364 | 18.225 | 1.29 ± 0.03 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr9 ../TE/HG00258.REF_chr9.bam` | 13.851 ± 0.166 | 13.699 | 14.168 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr9.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr9.bed` | 58.399 ± 1.429 | 56.511 | 60.452 | 4.22 ± 0.11 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr9.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr9.bed` | 121.349 ± 1.265 | 119.776 | 123.549 | 8.76 ± 0.14 |
