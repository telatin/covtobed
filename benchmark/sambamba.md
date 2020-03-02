# Benchmark including sambamba

*sambamba* can perform coverage statics calculation, but not producing a BED output. 

The "base" output will print a line for each position, while the "window" output will use a sliding window approach.

For the records, this benchmark has been run using the "example1.bam" target enrichment panel.

| Command | Mean [s] | Min [s] | Max [s] | Relative |
|:---|---:|---:|---:|---:|
| `mosdepth -F 4 -x -t 4 mos4_ example1.bam` | 76.689 ± 2.354 | 73.243 | 79.350 | 52.47 ± 1.75 |
| `mosdepth -F 4 -x      mos1_ example1.bam` | 77.538 ± 2.387 | 74.845 | 81.726 | 53.05 ± 1.78 |
| `covtobed example1.bam > covtobed_example1.bed` | 1.462 ± 0.019 | 1.436 | 1.489 | 1.00 |
| `sambamba depth base -F 'not unmapped' example1.bam > sambamba.txt.bed` | 4.376 ± 0.297 | 3.952 | 4.768 | 2.99 ± 0.21 |
| `sambamba depth window -F 'not unmapped' -w 1000 example1.bam > sambamba2.txt.bed` | 12.646 ± 0.624 | 11.966 | 13.425 | 8.65 ± 0.44 |
| `bedtools genomecov -bga -ibam example1.bam > bedtools_example1.bed` | 35.298 ± 3.425 | 32.689 | 40.519 | 24.15 ± 2.36 |
