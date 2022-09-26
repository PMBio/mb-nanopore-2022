#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../genome/hg38.fa
MM=/g/solexa/bin/genomesNew/mm10/mm10.fa
if [ -f ${HG} ]
then
    if [ -f ${MM} ]
    then
	if [ ! -f ref.fa ]
	then
	    cat ${HG} > ref.fa
	    cat ${MM} | sed 's/^>/>mouse/' >> ref.fa
	fi
	for F in *.fastq.gz
	do
	    if [ -f ${F} ]
	    then
		ID=`echo ${F} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
		if [ ! -f ${ID}.bam ]
		then
		    echo ${ID}
	    	    minimap2 --split-prefix tmpidx -t 6 -a -x map-ont -L ref.fa ${F} | samtools sort -o ${ID}.bam
		    samtools index ${ID}.bam
		fi
	    fi
	done
	rm ref.fa
    fi
fi

# Subset to human
for BAM in *.bam
do
    samtools view -b ${BAM} chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX > test.bam
    mv test.bam ${BAM}
    samtools index ${BAM}
done

