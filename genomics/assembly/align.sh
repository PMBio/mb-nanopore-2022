#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../genome/hg38.fa
ID=primary
if [ -f ${HG} ]
then
    if [ -f primary/assembly.fasta ]
    then
	# Against hg38
	minimap2 -t 6 -a -x asm5 -L ${HG} primary/assembly.fasta | samtools sort -o ${ID}.hg38.bam
	samtools index ${ID}.hg38.bam
    fi
fi
