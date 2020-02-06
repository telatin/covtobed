# Benchmark

## How to run

The benchmark has been performed using [hyperfine](https://github.com/sharkdp/hyperfine), installable via conda.
To download the datasets use the `get_datasets.sh` script. 
The script will download two target enrichment panels (BAM).
If invoked with a parameter, will also download an exome from the 1000 Genomes Project.

## Results

### Straming speed (not saving to disk)

_mosdepth_ is not included as it only saves to file. Other commands were redirected to `/dev/null`.

| Command | Mean [ms] | Min [ms] | Max [ms] | Relative |
|:---|---:|---:|---:|---:|
| `covtobed example1.bam` | 986.2 ± 9.5 | 977.1 | 1003.5 | 1.00 |
| `bedtools genomecov -bga -ibam example1.bam` | 33596.2 ± 1550.4 | 31982.4 | 35700.0 | 34.07 ± 1.61 |


### Saving the output to disk

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| ` mosdepth mosd2_ example2.bam` | 68.555 ± 0.765 | 67.672 | 68.999 | 33.57 ± 0.49 |
| `covtobed example2.bam > covt2_example2.bed` | 2.042 ± 0.019 | 2.024 | 2.061 | 1.00 |
| `bedtools genomecov -bga -ibam example2.bam > bedt2_example2.bed` | 35.986 ± 0.167 | 35.844 | 36.169 | 17.62 ± 0.18 |
