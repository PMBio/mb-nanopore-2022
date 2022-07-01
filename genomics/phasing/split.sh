#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate phasing

HG=${BASEDIR}/../genome/hg38.fa

for TYPE in rephase final
do
    mkdir split_${TYPE}
    bcftools view -s blood ${TYPE}.af.bcf | bcftools annotate -O b -o split_${TYPE}/blood.phased.bcf -x INFO,^FORMAT/GT -
    bcftools index split_${TYPE}/blood.phased.bcf

    # Split reads
    for BAM in ${BASEDIR}/../alignment/ont/*.bam
    do
	if [ -f ${BAM} ]
	then
	    ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	    echo ${ID}
	    if [ ${ID} == "Germline" ]; then continue; fi
	    samtools view -b ${BAM} `seq 1 22 | sed 's/^/chr/' | tr '\n' ' '` > tmp.bam
	    samtools index tmp.bam
	    /opt/dev/alfred/bin/alfred split -r ${HG} -s blood -v split_${TYPE}/blood.phased.bcf -p split_${TYPE}/${ID}.h1.bam -q split_${TYPE}/${ID}.h2.bam tmp.bam
	    rm tmp.bam tmp.bam.bai
	fi
    done
done
