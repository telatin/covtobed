| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_unmapped ../TE/HG00258.REF_unmapped.bam` | 1.387 ± 0.073 | 1.277 | 1.479 | 1.00 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_unmapped ../TE/HG00258.REF_unmapped.bam` | 1.430 ± 0.069 | 1.360 | 1.549 | 1.03 ± 0.07 |
| `covtobed ../TE/HG00258.REF_unmapped.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_unmapped.bed` | 9.131 ± 0.131 | 9.003 | 9.333 | 6.58 ± 0.36 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_unmapped.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_unmapped.bed` | 57.349 ± 0.506 | 56.851 | 58.155 | 41.35 ± 2.20 |
