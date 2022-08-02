#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../../genome/hg38.fa

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ! -f ${ID}.json.gz ]
	then
	    echo ${ID}
	    alfred qc -r ${HG} -j ${ID}.json.gz -o ${ID}.tsv.gz ${BAM}
	fi
    fi
done
