#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate visualization

HG=${BASEDIR}/../../genome/hg38.fa
if [ -f ${HG} ]
then
    cd ${BASEDIR}
    for ALIGN in ../alignment/ont/*.bam
    do
	if [ -f ${ALIGN} ]
	then
	    ID=`echo ${ALIGN} | sed 's/^.*\///' | sed 's/H021-//' | sed 's/-.*$//'`
	    echo ${ID}
	    if [ ${ID} == "D7K5QE" ]
	    then
		samtools view ${ALIGN} chr12:103700218-103700284 | cut -f 1 | sort | uniq -d > reads
	    else
		samtools view ${ALIGN} chr16:63892950-63893800 chr16:57105500-57105800 | cut -f 1 | sort | uniq -d > reads
	    fi
	    wally dotplot -g ${HG} -R reads -s 13000 ${ALIGN}
	    rm reads
	fi
    done
fi
