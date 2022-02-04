#!/bin/bash

for CHR in 3
do
    for METHOD in whatshap shapeit
    do
	INFILE=../${METHOD}.af.bcf
	if [ -f ${INFILE} ]
	then
	    #[1]POS,[2]tumor:GT,[3]tumor:AD,[4]relapse:GT,[5]relapse:AD,[6]blood:GT,[7]blood:AD
	    bcftools view ${INFILE} chr${CHR} | bcftools query -f "%POS[\t%GT\t%AD]\n" - | awk '$6=="0|1" || $6=="1|0"' | sed 's/,/\t/g' | awk '$3+$4>0 {print $1"\t"$4/($3+$4)"\t"$8;}' | grep "0|1" | awk '{print $0"\t"(1-$2)"\t"$2;}' > partA
	    bcftools view ${INFILE} chr${CHR} | bcftools query -f "%POS[\t%GT\t%AD]\n" - | awk '$6=="0|1" || $6=="1|0"' | sed 's/,/\t/g' | awk '$3+$4>0 {print $1"\t"$4/($3+$4)"\t"$8;}' | grep "1|0" | awk '{print $0"\t"$2"\t"(1-$2);}' > partB
	    sort -k1,1n partA partB > chr${CHR}.${METHOD}
	    Rscript plot.R chr${CHR}.${METHOD}
	    rm partA partB chr${CHR}.${METHOD}
	fi
    done
done
