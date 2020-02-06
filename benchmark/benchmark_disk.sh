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

	if [ ! -z ${1+x} ];
	then
		# skip exome
		if [ "$TAG" == "HG00258" ]; then
			echo "Skipping exome: specify any parameter to include it";
			exit;
		fi
	fi

	if [ ! -e "$FILE" ]; then
	    echo "$FILE not found. Try re-downloading the datasets"
	    exit 2
	fi

	echo "$TAG [$FILE]"
	hyperfine --warmup 1 --min-runs 3 --prepare "rm $FILE.bai mosd_* *.bed || true" --cleanup 'rm *.XXX || true' \
		--export-csv benchmark_$TAG.csv --export-markdown benchmark_$TAG.md \
		"samtools index $FILE && mosdepth \"$DIR\"/mosd_$TAG $FILE && gunzip \"$DIR\"/mosd_$TAG.per-base.bed.gz" \
		"covtobed $FILE > \"$DIR\"/covt_$TAG.bed" \
		"bedtools genomecov -bga -ibam $FILE > \"$DIR\"/bedt_$TAG.bed" 

done
