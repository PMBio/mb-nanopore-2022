#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate phasing

FILES=""
for CHR in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do
    bcftools view -O b -o chr${CHR}.rephase.bcf shapeit.af.bcf chr${CHR} 
    for SNPS in `seq 20 40 100 | sort -nr`
    do
	for T in tumor relapse
	do
	    for I in `seq 1 5`
	    do
		echo "chr${CHR}, iteration..." ${I} ${SNPS} ${T}
		mv chr${CHR}.rephase.bcf chr${CHR}.rephase.in.bcf
		python phase.py -s ${SNPS} -v chr${CHR}.rephase.in.bcf -g blood -t ${T} -o chr${CHR}.rephase.bcf
		rm chr${CHR}.rephase.in.bcf
	    done
	done
    done
    bcftools index chr${CHR}.rephase.bcf
    FILES=${FILES}" chr${CHR}.rephase.bcf"
done

# Concatenate
bcftools concat -O b -o rephase.af.bcf ${FILES}
bcftools index rephase.af.bcf
rm chr*rephase.bcf*
