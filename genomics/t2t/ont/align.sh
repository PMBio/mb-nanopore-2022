#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate align

HG=${BASEDIR}/../../genome/t2t.fa
if [ -f ${HG} ]
then
    cd ${BASEDIR}
    for F in ${BASEDIR}/../../fastq/ont/*.fastq.gz
    do
	if [ -f ${F} ]
	then
	    ID=`echo ${F} | sed 's/^.*\///' | sed 's/.fastq.gz$//'`
	    if [ ! -f ${ID}.bam ]
	    then
		echo ${ID}
	    	minimap2 -t 24 -a -x map-ont -L ${HG} ${F} | samtools sort -@ 12 -m 4608M -o ${ID}.bam
		samtools index -@ 12 ${ID}.bam
	    fi
	fi
    done
fi

# SV calling
/opt/dev/delly/src/delly lr -y ont -g ../../genome/t2t.fa -o delly.bcf *.bam

delly cnv -a -g ${HG} -m ${MAP} -c ${ID}.cov.gz -o ${ID}.bcf ${BAM}
