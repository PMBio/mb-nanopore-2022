#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Find telomeric motifs
for BAM in ${BASEDIR}/../../alignment/illumina/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	/opt/dev/alfred/src/alfred telmotif -r ${HG} -o ${ID}.bed ${BAM}
    fi
done

# Filter
for TYPE in tumor relapse
do
    echo ${TYPE}
    cat ${TYPE}.bed  | grep -P "^chr[0-9X]*\t" | awk '$4>=5' > ${TYPE}.selected
    bedtools intersect -a ${TYPE}.selected -b blood.bed -wao | awk '$5=="."' | cut -f 1-4
    rm ${TYPE}.selected
done
