#!/bin/bash

# This scripts allows to catch alterations of the code producing a different output from expected results
# The typical test is a comparison of number of lines produced
# Runs automatically with TravisCI, and is executed from the root of the repository.
# REQUIRES: ./test/demo.bam

set -exo

# Compilation success, checking that --version emits the expected "progname" string
echo -n "Compiled binary prints version: "
if [ $(./covtobed --version  | grep covtobed | wc -l ) -eq "1" ];
then
        echo PASS
fi

# Testing that -m MIN produces an output, and it fits the expectation for demo.bam (n. lines)
echo -n "Minimum coverage, expected BED lines check: "
if [ $(./covtobed -m 15 test/demo.bam  | wc -l) -eq "12" ];
then
	echo PASS
fi

# Checking thath --physical-coverage will work, and it fits the expected number of lines
echo -n "Physical coverage, expected BED lines check: "
if [ $(./covtobed --physical-coverage test/mp.bam  | wc -l) -eq "136" ];
then
	echo PASS
fi

# Checking stranded output: it should produce content in the fifth column of the bed file
echo -n "Stranded output, testing column #5: "
if [ $(./covtobed --out test/demo.bam | cut -f 5 | sort -u | wc -l) -eq "10" ];
then
	echo PASS
fi

# Checking the "counts" output (counting the lines containing a ">")
echo -n "Testing 'counts' format (printed headers): "
if [ $(./covtobed --format counts test/demo.bam | grep '>' | wc -l) -eq "2" ];
then
	echo PASS
fi
echo -n "Testing 'counts' format (printed lines): "
if [ $(./covtobed --format counts test/demo.bam | grep -v '>' | wc -l) -eq "202" ];
then
        echo PASS
fi

# Checking BED output with reference output file
echo -n "Checking identity of BED output with pre-calculated: "
./covtobed test/demo.bam > test/output.test
if [ $(diff test/output.test test/output.bed | wc -l) -eq "0" ];
then
        echo PASS
	rm test/output.test
fi
