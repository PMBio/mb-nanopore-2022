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
	samtools calmd -b ${BAM} ${HG} > tmp.bam 2> /dev/null
	samtools index tmp.bam
	sniffles -m tmp.bam -v ${ID}.vcf
	rm tmp.bam*
    fi
done
