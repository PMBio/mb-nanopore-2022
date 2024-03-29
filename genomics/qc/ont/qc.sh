#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../genome/hg38.fa

for FQ in ${BASEDIR}/../../fastq/ont/*.fastq.gz
do
    if [ -f ${FQ} ]
    then
	ID=`echo ${FQ} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
	echo ${ID}
	NanoStat --fastq ${FQ} > ${ID}.nanostat
    fi
done

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	alfred qc -r ${HG} -j ${ID}.json.gz -o ${ID}.tsv.gz ${BAM}
    fi
done
