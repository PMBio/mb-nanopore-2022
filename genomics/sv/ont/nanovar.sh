#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then

	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	nanovar -x ont -t 8 -f hg38 ${BAM} ${HG} ${ID}_dir
    fi
done
