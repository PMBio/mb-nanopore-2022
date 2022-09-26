#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../genome/hg38.fa
MAP=Homo_sapiens.GRCh38.dna.primary_assembly.fa.e1.r101.blacklist.gz
WIN=500000

for BAM in *.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	echo ${ID}
	delly cnv -i ${WIN} -j ${WIN} -w ${WIN} -a -g ${HG} -m ${MAP} -c ${ID}.cov.gz -o ${ID}.bcf ${BAM}
    fi
done
