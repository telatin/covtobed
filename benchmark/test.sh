#!/bin/bash

set -euxo pipefail

# Uses: 'hyperfine', available with 'conda install -y  -c conda-forge hyperfine'

for FILE in ../pannelli/16.bam ../pannelli/17.bam ../exome/HG00258.bam;
do
	TAG=$(basename $FILE | cut -f1 -d.)
	hyperfine --warmup 4 --min-runs 15 --cleanup 'rm *.bed || true' \
		--export-csv benchmark_$TAG.csv --export-markdown benchmark_$TAG.md \
		"covtobed $FILE > /dev/null" \
		"bedtools genomecov -bga -ibam $FILE > /dev/null"
done
