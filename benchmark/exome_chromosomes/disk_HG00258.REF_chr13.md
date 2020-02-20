| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr13 ../TE/HG00258.REF_chr13.bam` | 13.006 ± 0.390 | 12.620 | 13.740 | 1.28 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr13 ../TE/HG00258.REF_chr13.bam` | 10.171 ± 0.238 | 9.852 | 10.421 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr13.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr13.bed` | 40.300 ± 1.262 | 39.388 | 42.717 | 3.96 ± 0.15 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr13.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr13.bed` | 94.534 ± 0.800 | 93.509 | 95.422 | 9.29 ± 0.23 |
