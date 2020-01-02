set -e
echo -n "Compiled: "
if [ $(./covtobed --version  | grep covtobed | wc -l ) -eq 1 ];
then
        echo PASS
fi

echo -n "Minimum coverage: "
if [ $(./covtobed -m 15 test/demo.bam  | wc -l) -eq "12" ];
then
	echo PASS
fi

echo -n "Physical coverage: "
if [ $(./covtobed --physical-coverage test/mp.bam  | wc -l) -eq 136 ];
then
	echo PASS
fi

