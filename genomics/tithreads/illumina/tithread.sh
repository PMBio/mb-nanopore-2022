#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Discovery templated insertions
for BAM in ${BASEDIR}/../../alignment/illumina/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ${ID} == "blood" ]; then continue; fi
	/opt/dev/rayas/src/rayas call -g ${HG} -o ${ID}.bed -m ${BASEDIR}/../../alignment/illumina/blood.bam ${BAM}
    fi
done

# Filter for self-concatenation and high copy-number
echo -e "sample\tinstances\tloci" > stats.tsv
for TYPE in tumor relapse
do
    INST=`tail -n +2 ${TYPE}.bed | awk '$5>=3 && $6>=50 && $7>=10' | cut -f 8 | sort | uniq  | wc -l | cut -f 1`
    NUM=0
    if [ ${INST} -gt 0 ]
    then
	for CLUSTID in `tail -n +2 ${TYPE}.bed | awk '$5>=3 && $6>=50 && $7>=10' | cut -f 8 | sort | uniq`
	do
	    SEGMENTS=`cat ${TYPE}.bed | awk '$8=="'${CLUSTID}'"' | cut -f 7 | awk '{SUM+=$1; BASE+=2;} END {print (SUM-BASE);}'`
	    echo ${SEGMENTS} | sed "s/^/${TYPE}\t${CLUSTID}\t/" >> stats.ti.seglen
	    cat ${TYPE}.bed | awk '$8=="'${CLUSTID}'"' | sed "s/^/${TYPE}\t/" >> stats.ti.tsv
	    NUMCLUST=`cat ${TYPE}.bed | awk '$8=="'${CLUSTID}'"' | cut -f 1-3 | sort | uniq  | wc -l | cut -f 1`
	    NUM=`expr ${NUM} + ${NUMCLUST}`
	done
    fi		
    echo -e "${TYPE}\t${INST}\t${NUM}" >> stats.tsv
done
