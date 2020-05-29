# Example output

Output of *mosdepth*, *covtobed* and *bedtools genomecov* from the ../demo.bam file

```bash
covtobed ../demo.bam > covtobed.bed
bedtools genomecov -bga  -ibam ../demo.bam > bedtools.bed
mosdepth mos ../demo.bam && zcat mos.per-base.bed.gz > mosdepth.bed && rm mos.*
```

For the given example, *covtobed* and *mosdepth* produce an identical output, while *bedtools* will print the first line of covtobed at the end of its output. 

See [Wiki: implementation](https://github.com/telatin/covtobed/wiki/covtobed) to understand the potential source of difference between the tools.
