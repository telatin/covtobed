| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrEBV ../TE/HG00258.REF_chrEBV.bam` | 1.571 ± 0.067 | 1.498 | 1.693 | 1.03 ± 0.05 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_HG00258.REF_chrEBV ../TE/HG00258.REF_chrEBV.bam` | 1.524 ± 0.040 | 1.478 | 1.588 | 1.00 |
| `covtobed ../TE/HG00258.REF_chrEBV.bam > "/local/giovanni/covtobed/benchmark"/covt2_HG00258.REF_chrEBV.bed` | 1.669 ± 0.062 | 1.584 | 1.777 | 1.10 ± 0.05 |
| `bedtools genomecov -bga -ibam ../TE/HG00258.REF_chrEBV.bam > "/local/giovanni/covtobed/benchmark"/bedt2_HG00258.REF_chrEBV.bed` | 39.349 ± 1.039 | 38.270 | 40.950 | 25.83 ± 0.96 |
