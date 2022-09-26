#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../genome/hg38.fa
MM=/g/solexa/bin/genomesNew/mm10/mm10.fa

for FQ in *.fastq.gz
do
    if [ -f ${FQ} ]
    then
	ID=`echo ${FQ} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
	NanoStat --fastq ${FQ} > ${ID}.nanostat
    fi
done

if [ ! -f ref.fa ]
then
    cat ${HG} > ref.fa
    cat ${MM} | sed 's/^>/>mouse/' >> ref.fa
    samtools faidx ref.fa
fi

for BAM in *.bam
do
    if [ -f ${BAM} ]
    then
	ID=`echo ${BAM} | sed 's/^.*\///' | sed 's/.bam$//'`
	alfred qc -r ref.fa -j ${ID}.json.gz -o ${ID}.tsv.gz ${BAM}
    fi
done
rm ref.fa*
