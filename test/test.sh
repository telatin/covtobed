#!/bin/bash

# This scripts allows to catch alterations of the code producing a different output from expected results
# The typical test is a comparison of number of lines produced
# Runs automatically with TravisCI, and is executed from the root of the repository.
# REQUIRES: ./test/demo.bam
#           ./test/mock.bam


if [ ! -e "test/demo.bam" ]; then
	echo "WARNING: running this test from a bad location"
	DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
	cd "$DIR"/..
fi

if [ ! -e "./covtobed" ]; then
	echo "WARNING:"
	echo "Binary not found: tring to use pre-compiled...";
	if [ `uname` == 'Darwin' ]; then
		echo " - Copying macOS binary"  || echo " ERROR: pre-compiled binary not working"
		cp ./binaries/covtobed_mac ./covtobed
		covtobed --version
	else
		echo " - Trying Ubuntu binary"
		cp ./binaries/covtobed ./covtobed
		covtobed --version || echo " ERROR: pre-compiled binary not working"
	fi
fi

./covtobed --version

set -eou pipefail
# Compilation success, checking that --version emits the expected "progname" string
echo -n " - Compiled binary prints version: "
if [ $(./covtobed --version  | grep covtobed | wc -l ) -eq "1" ];
then
        echo PASS 1
fi

# Testing that -m MIN produces an output, and it fits the expectation for demo.bam (n. lines)
echo -n " - Minimum coverage, expected BED lines check: "
if [ $(./covtobed -m 15 test/demo.bam  | wc -l) -eq "12" ];
then
	echo PASS 2
else
	echo FAIL
	exit 1
fi

# Checking thath --physical-coverage will work, and it fits the expected number of lines
echo -n " - Physical coverage, expected BED lines check: "
if [ $(./covtobed --physical-coverage test/mp.bam  | wc -l) -eq "136" ];
then
	echo PASS 3
else
	echo FAIL
	exit 1
fi

# Checking stranded output: it should produce content in the fifth column of the bed file
echo -n " - Stranded output, testing column #5: "
if [ $(./covtobed --out test/demo.bam | cut -f 5 | sort -u | wc -l) -eq "10" ];
then
	echo PASS 4
else
	echo FAIL
	exit 1
fi

# Checking the "counts" output (counting the lines containing a ">")
echo -n " - Testing 'counts' format (printed headers): "
if [ $(./covtobed --format counts test/demo.bam | grep '>' | wc -l) -eq "2" ];
then
	echo PASS 5
else
	echo FAIL
	exit 1
fi
echo -n " - Testing 'counts' format (printed lines): "
if [ $(./covtobed --format counts test/demo.bam | grep -v '>' | wc -l) -eq "202" ];
then
        echo PASS 6
else
	echo FAIL
	exit 1
fi

# Checking BED output with reference output file
echo -n " - Checking identity of BED output with pre-calculated: "
./covtobed test/demo.bam > test/output.test
if [ $(diff test/output.test test/output.bed | wc -l) -eq "0" ];
then
        echo PASS 7
	rm test/output.test
else
	echo FAIL
	exit 1
fi

## Synthetic BAM test
echo -n " - Checking computed coverage for a synthetic BAM file: "
./covtobed -m 1 test/mock.bam > test/output.test
if [ $(diff test/output.test test/mock.bed | wc -l) -eq "0" ];
then
        echo PASS 8
	rm test/output.test
else
	echo FAIL
	exit 1
fi

## Filter non valid alignments
echo -n " - Checking filtering of invalid alignments: "
if [ $(./covtobed -m 1 -d test/filtered.bam | wc -l) -eq "2" ] ; then
	echo -n "PASS 9,"
else
	echo FAIL
	exit 1
fi
if [ $(./covtobed -m 1 test/filtered.bam | wc -l) -eq "5" ] ; then
	echo 10
else
	echo "FAIL: $(./covtobed -m 1 -d test/filtered.bam | wc -l)"
	exit 1
fi

# A synthetic BAM containing reference name {n}X that should only print one reagion 
# covered exactly {n}X times
echo  " - Checking artificial coverage values:"
./covtobed test/test_cov.bam -m 1 |cut -f 1,4| while read LINE;
do
	echo $LINE | perl -ne '($exp, $cov)=split /\s+/, $_; if ("$exp" ne "${cov}X") {
		die "$exp != $cov\n";
	} else {
		print "\t#OK expecting $exp, ${cov}X found\n";
	}'
done
echo "ALL TESTS: PASSED"
