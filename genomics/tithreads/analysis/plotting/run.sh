#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate assembly

# Plot individual reads, as self-alignments and against the templated insertion source loci
for READ in 0645e9e8-3cc3-4f9f-b6f1-382412b524b6 5cce2684-dd22-40c5-8e2b-97164edfb753 b6678d21-0115-4a5a-81b8-eb5e8e269890
do
    zcat ${BASEDIR}/../../ont/Primary_tumor.fa.gz | grep -w "${READ}" -A 1 > ${READ}.fa
    ./selfalign.sh ${READ}.fa

    # Plot against the 2 templated insertion threads
    zcat ${BASEDIR}/../../ont/Primary_tumor.match.gz | grep -w "${READ}" > matches.${READ}.tsv
    for REGIONS in ${BASEDIR}/../overlap_ont_illumina/core*.bed
    do
	CLUSTID=`echo ${REGIONS} | sed 's/^.*\///' | sed 's/^core.//' | sed 's/.bed//'`
	bedtools intersect -a matches.${READ}.tsv -b ${REGIONS} -wao | awk '$9!="."' | sort -k1,1V -k2,2n | uniq > matches.${READ}.${CLUSTID}
	if [ `cat matches.${READ}.${CLUSTID} | wc -l` -gt 0 ]
	then
	    Rscript projection.R matches.${READ}.${CLUSTID}
	fi
	rm matches.${READ}.${CLUSTID}
    done
    rm matches.${READ}.tsv ${READ}.fa
done
