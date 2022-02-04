#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

for BAM in ${BASEDIR}/../../alignment/illumina/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	/opt/dev/alfred/src/alfred telmotif -r ${HG} -o ${ID}.bed ${BAM}
    fi
done
