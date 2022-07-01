#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa
MAP=Homo_sapiens.GRCh38.dna.primary_assembly.fa.e1.r101.blacklist.gz

for BAM in ${BASEDIR}/../../alignment/illumina/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	delly cnv -a -g ${HG} -m ${MAP} -c ${ID}.cov.gz -o ${ID}.bcf ${BAM}
    fi
done
