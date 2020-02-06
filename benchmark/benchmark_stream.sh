#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

set -euo pipefail

# Uses: 'hyperfine', available with 'conda install -y  -c conda-forge hyperfine'

command -v hyperfine || exit 1

if [ ! -e "$DIR"/../TE/example1.bam ]; then
    echo "Download the input datasets first"
    exit
fi
for FILE in "$DIR"/../TE/example1.bam "$DIR"/../TE/example2.bam "$DIR"/../TE/HG00258.bam;
do
	TAG=$(basename $FILE | cut -f1 -d.)

	if [ ! -e "$FILE" ]; then
	    echo "$FILE not found. Try re-downloading the datasets"
	    exit 2
	fi
	if [ ! -z ${1+x} ];
	then
		# skip exome
		if [ "$TAG" == "HG00258" ]; then
			echo "Skipping exome: specify any parameter to include it";
			exit;
		fi
	fi
	echo "$TAG [$FILE]"
	hyperfine --warmup 1 --min-runs 6 --cleanup 'rm *.bed || true' \
		--export-csv stream/benchmark_$TAG.csv --export-markdown stream/benchmark_$TAG.md \
		"samtools depth -a $FILE > /dev/null" \
		"covtobed $FILE > /dev/null" \
		"bedtools genomecov -bga -ibam $FILE > /dev/null"
done
