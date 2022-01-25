#!/bin/bash

if [ $# -lt 1 ]
then
    echo "Usage: $0 chr1:1000-2000"
    exit -1
fi

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate variants

for BAM in ${BASEDIR}/*/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	samtools view -b ${BAM} $@ | samtools sort -o ${ID}.bam -
	samtools index ${ID}.bam
    fi
done
