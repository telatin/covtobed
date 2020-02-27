
### Human panel - Straming speed (not saving to disk)

_mosdepth_ is not included as it only saves to file. Other commands were redirected to `/dev/null`.
See also [example2.bam benchmark](stream/benchmarkStream_example2.md). 

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `bedtools genomecov -bga -ibam example1.bam` | 34.717 ± 0.952 | 33.325 | 36.090 | 33.17 ± 1.79 |
| `covtobed example1.bam` | 1.046 ± 0.049 | 0.982 | 1.096 | 1.00 |
