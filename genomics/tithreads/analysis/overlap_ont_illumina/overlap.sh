#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate variants

if [ -f ${BASEDIR}/../../ont/stats.ti.tsv ]
then
    for CLUSTID in `cat ../../ont/stats.ti.tsv  | cut -f 9 | sort | uniq`
    do
	echo "Primary tumor, cluster ${CLUSTID}"
	cat ${BASEDIR}/../../ont/stats.ti.tsv | awk '$9=="'${CLUSTID}'"' | cut -f 2-4,8 | sort -k1,1V -k2,2n | uniq > ont.bed
	cat ${BASEDIR}/../../illumina/stats.ti.tsv | cut -f 2-4,8 | sort -k1,1V -k2,2n | uniq > illumina.bed
	bedtools  intersect -a ont.bed -b illumina.bed -wao
	rm ont.bed illumina.bed 
    done
fi
