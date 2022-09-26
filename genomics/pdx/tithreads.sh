#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../genome/hg38.fa

for BAM in *.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ${ID} == "BT084" ]
	then
	    CTRL=Prim
	else
	    CTRL=BT084
	fi
	lorax tithreads -o ${ID} -g ${HG} -m ${CTRL}.bam ${BAM}
    fi
done
