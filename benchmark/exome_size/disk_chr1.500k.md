| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -x "/local/giovanni/covtobed/benchmark"/mosd2_chr1.500k ../TE/chr1.500k.bam` | 8.659 ± 0.265 | 8.386 | 9.148 | 7.68 ± 0.27 |
| `mosdepth -x -t 4 "/local/giovanni/covtobed/benchmark"/mosd2_chr1.500k ../TE/chr1.500k.bam` | 8.337 ± 0.347 | 7.943 | 8.908 | 7.39 ± 0.34 |
| `covtobed ../TE/chr1.500k.bam > "/local/giovanni/covtobed/benchmark"/covt2_chr1.500k.bed` | 1.127 ± 0.021 | 1.100 | 1.154 | 1.00 |
| `bedtools genomecov -bga -ibam ../TE/chr1.500k.bam > "/local/giovanni/covtobed/benchmark"/bedt2_chr1.500k.bed` | 35.909 ± 1.480 | 34.224 | 38.168 | 31.85 ± 1.44 |
