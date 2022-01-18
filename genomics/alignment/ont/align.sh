#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
if [ -f ${HG} ]
then
    cd ${BASEDIR}
    for F in ${BASEDIR}/../../fastq/ont/*.fastq.gz
    do
	if [ -f ${F} ]
	then
	    ID=`echo ${F} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
	    echo ${ID}
	    minimap2 -t 6 -a -L -x map-ont -L ${HG} ${F} | samtools sort -o ${ID}.bam
	    samtools index ${ID}.bam
	fi
    done
fi
