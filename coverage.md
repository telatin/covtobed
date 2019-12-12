# coverage (covtobed legacy ancestor)

**coverage** was a program developed by Giovanni Birolo, that we included in the [Docker image](https://hub.docker.com/r/andreatelatin/covtobed) of **covtobed**.

## Synopsis

```
Usage: coverage [options]

Computes coverage from alignments

Options:
  -h, --help            show this help message and exit
  -i BAM, --input=BAM   a sorted bam file with the input alignments
  -t BED, --target=BED  skip regions outside the intervals in the sorted bed
file BED
  --no-out              do not output coverage
  --wig                 outputs coverage in wig format
  -b, --bed             outputs coverage in bed format, one line per base
  -r, --regions         outputs intervals where coverage is in the range given
                        by --min-cov and --max-cov, in bed format
  -f, --filtered        outputs coverage only for bases read from standard input
  -p, --physical-cov    computes physical coverage for mate pairs
  --min-qual=MIN        ignore alignments whose quality is less than MIN
(default: 0)
  -m MIN, --min-cov=MIN 
                        skip regions where coverage is less than MIN (default: 1)
  -M MAX, --max-cov=MAX 
                        skip regions where coverage is more than MAX (default: +INF)
  -s, --summary         print stats to standard error
  --hide-summary        do not print stats to standard error
  -q, --quiet           do not print progress
  --show-progress       print progress (the default if standard error is a terminal
  --hide-progress       do not print progress (the default if standard error is
                        not a terminal)
```
