| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr6 ../TE/HG00258.REF_chr6.bam` | 22.374 ± 0.301 | 22.039 | 22.928 | 1.32 ± 0.04 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr6 ../TE/HG00258.REF_chr6.bam` | 16.896 ± 0.494 | 16.305 | 17.702 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr6.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr6.bed` | 78.617 ± 2.604 | 75.033 | 82.593 | 4.65 ± 0.21 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr6.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr6.bed` | 147.321 ± 3.141 | 142.801 | 150.217 | 8.72 ± 0.32 |
