| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrY ../TE/HG00258.REF_chrY.bam` | 3.444 ± 0.102 | 3.322 | 3.570 | 1.53 ± 0.06 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrY ../TE/HG00258.REF_chrY.bam` | 3.165 ± 0.083 | 3.071 | 3.272 | 1.40 ± 0.05 |
| `covtobed ../TE/HG00258.REF_chrY.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chrY.bed` | 2.257 ± 0.057 | 2.205 | 2.370 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chrY.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chrY.bed` | 39.259 ± 0.417 | 38.739 | 39.800 | 17.39 ± 0.48 |
