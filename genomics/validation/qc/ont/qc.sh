#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../../genome/hg38.fa

for FQ in ${BASEDIR}/../../fastq/ont/*.fastq.gz
do
    if [ -f ${FQ} ]
    then
	ID=`echo ${FQ} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
	if [ ! -f ${ID}.nanostat ]
	then
	    echo ${ID}
	    NanoStat --fastq ${FQ} > ${ID}.nanostat
	fi
    fi
done

for BAM in ${BASEDIR}/../../alignment/ont/*.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	if [ ! -f ${ID}.json.gz ]
	then
	    echo ${ID}
	    alfred qc -r ${HG} -j ${ID}.json.gz -o ${ID}.tsv.gz ${BAM}
	fi
    fi
done
