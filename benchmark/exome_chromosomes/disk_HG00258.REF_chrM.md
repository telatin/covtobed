| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrM ../TE/HG00258.REF_chrM.bam` | 1.681 ± 0.045 | 1.623 | 1.724 | 1.13 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrM ../TE/HG00258.REF_chrM.bam` | 1.543 ± 0.072 | 1.470 | 1.654 | 1.04 ± 0.06 |
| `covtobed ../TE/HG00258.REF_chrM.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chrM.bed` | 1.483 ± 0.049 | 1.441 | 1.573 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chrM.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chrM.bed` | 39.423 ± 1.034 | 38.305 | 40.769 | 26.58 ± 1.12 |
