| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr1 ../TE/HG00258.REF_chr1.bam` | 35.681 ± 0.650 | 35.154 | 36.974 | 1.41 ± 0.03 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr1 ../TE/HG00258.REF_chr1.bam` | 25.347 ± 0.362 | 25.044 | 25.940 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr1.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr1.bed` | 128.576 ± 3.142 | 122.704 | 131.607 | 5.07 ± 0.14 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr1.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr1.bed` | 223.619 ± 6.603 | 217.182 | 235.391 | 8.82 ± 0.29 |
