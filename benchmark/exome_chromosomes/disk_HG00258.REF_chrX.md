| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrX ../TE/HG00258.REF_chrX.bam` | 20.624 ± 0.675 | 19.881 | 21.792 | 1.35 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrX ../TE/HG00258.REF_chrX.bam` | 15.313 ± 0.301 | 14.888 | 15.653 | 1.00 |
| `covtobed ../TE/HG00258.REF_chrX.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chrX.bed` | 68.434 ± 2.198 | 64.544 | 70.533 | 4.47 ± 0.17 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chrX.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chrX.bed` | 135.661 ± 2.760 | 132.815 | 139.865 | 8.86 ± 0.25 |
