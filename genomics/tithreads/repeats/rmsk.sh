#!/bin/bash

REGIONS=tithreads

# Get repeat masker track
for HG in hg19 hg38
do
    echo ${HG}
    wget http://hgdownload.cse.ucsc.edu/goldenpath/${HG}/database/rmsk.txt.gz
    zcat rmsk.txt.gz | awk -v OFS="\t" '{ print $6, $7, $8, $12, $11, $10 }' | sort -k1,1V -k2,2n | uniq | gzip -c > rmsk.${HG}.bed.gz
    rm rmsk.txt.gz

    # Split by repeat class
    for REPEAT in `zcat rmsk.${HG}.bed.gz | cut -f 4 | grep -v "?" | sort | uniq`
    do
	echo ${REPEAT}
	zcat rmsk.${HG}.bed.gz | grep -w "${REPEAT}" | grep -v "^[^\t]*_" | gzip -c > rmsk.${REPEAT}.${HG}.bed.gz
	if [ `zcat rmsk.${HG}.bed.gz | head | wc -l` -gt 0 ]
	then
	    ./shuffle.sh ${HG} 1000 ${REGIONS}.${HG}.bed rmsk.${REPEAT}.${HG}.bed.gz
	fi
    done
done
