# Example output

Output of mosdepth, covtobed and bedtools genomecov from the ../demo.bam file


```
covtobed ../demo.bam > covtobed.bed
bedtools genomecov -bga  -ibam ../demo.bam > bedtools.bed
mosdepth mos ../demo.bam && zcat mos.per-base.bed.gz > mosdepth.bed && rm mos.*
```
