#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=../../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
for F in ../../fastq/ont/*.fastq.gz
do
    if [ -f ${F} ]
    then
	ID=`echo ${F} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
	echo ${ID}
	minimap2 -t 6 -a -L -x map-ont -L ${HG} ${F} | samtools sort -o ${ID}.minimap2.srt.bam
	samtools index ${ID}.minimap2.srt.bam
    fi
done
