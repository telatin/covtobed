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
	hyperfine --warmup 1 --min-runs 3 --prepare "rm  mosd_* *.bed || true" --cleanup 'rm *.XXX || true' \
		--export-csv benchmark2_$TAG.csv --export-markdown benchmark2_$TAG.md \
		" mosdepth \"$DIR\"/mosd2_$TAG $FILE" \
		"covtobed $FILE > \"$DIR\"/covt2_$TAG.bed" \
		"bedtools genomecov -bga -ibam $FILE > \"$DIR\"/bedt2_$TAG.bed" 

done
