| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr10 ../TE/HG00258.REF_chr10.bam` | 18.382 ± 0.328 | 18.007 | 18.905 | 1.34 ± 0.03 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chr10 ../TE/HG00258.REF_chr10.bam` | 13.698 ± 0.218 | 13.444 | 14.099 | 1.00 |
| `covtobed ../TE/HG00258.REF_chr10.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chr10.bed` | 64.620 ± 1.886 | 61.831 | 67.343 | 4.72 ± 0.16 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chr10.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chr10.bed` | 126.262 ± 2.543 | 121.461 | 128.303 | 9.22 ± 0.24 |
