# Benchmark

## How to run

The benchmark has been performed using [hyperfine](https://github.com/sharkdp/hyperfine), installable via conda.
To download the datasets use the `get_datasets.sh` script. 
The script will download two target enrichment panels (BAM).
If invoked with a parameter, will also download an exome from the 1000 Genomes Project.

## Results: Linux VM (64 Gb RAM, 4 cores)

### Straming speed (not saving to disk)

_mosdepth_ is not included as it only saves to file. Other commands were redirected to `/dev/null`.
See also [example1.bam benchmark](stream/benchmarkStream_example1.md). 

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed /local/giovanni/covtobed/benchmark/../TE/example2.bam > /dev/null` | 1.951 ± 0.080 | 1.872 | 2.074 | 1.00 |
| `bedtools genomecov -bga -ibam /local/giovanni/covtobed/benchmark/../TE/example2.bam > /dev/null` | 36.709 ± 0.878 | 35.784 | 38.083 | 18.82 ± 0.90 |



### Saving the output to disk

This is the test done saving to file. 
Note that _mosdepth_ will save the file compressed and indexed, thus requiring more time, 
and it's the only program tested supporting multithreading. 

See also [example1.bam benchmark](disk/benchmark2_example1.md).

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| ` mosdepth mosd2_ example2.bam` | 68.555 ± 0.765 | 67.672 | 68.999 | 33.57 ± 0.49 |
| `covtobed example2.bam` | 2.042 ± 0.019 | 2.024 | 2.061 | 1.00 |
| `bedtools genomecov -bga -ibam example2.bam` | 35.986 ± 0.167 | 35.844 | 36.169 | 17.62 ± 0.18 |

## Results: macOS 10.15, MacBook Pro (16-inch, 2019)

_mosdepth_ is not available for macOS. 
Streaming test.

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed /Users/telatina/git/covtobed/benchmark/../TE/example2.bam > /dev/null` | 1.314 ± 0.024 | 1.289 | 1.346 | 1.00 |
| `bedtools genomecov -bga -ibam /Users/telatina/git/covtobed/benchmark/../TE/example2.bam > /dev/null` | 18.190 ± 0.134 | 18.071 | 18.398 | 13.84 ± 0.27 |
