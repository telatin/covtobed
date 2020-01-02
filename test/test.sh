set -e

echo -n "Compiled binary prints version: "
if [ $(./covtobed --version  | grep covtobed | wc -l ) -eq 1 ];
then
        echo PASS
fi

echo -n "Minimum coverage, expected BED lines check: "
if [ $(./covtobed -m 15 test/demo.bam  | wc -l) -eq "12" ];
then
	echo PASS
fi

echo -n "Physical coverage, expected BED lines check: "
if [ $(./covtobed --physical-coverage test/mp.bam  | wc -l) -eq 136 ];
then
	echo PASS
fi
echo -n "Stranded output, testing column #5: "
if [ $(./covtobed --out test/demo.bam | cut -f 5 | sort -u | wc -l) -eq 10 ];
then
	echo PASS
fi
