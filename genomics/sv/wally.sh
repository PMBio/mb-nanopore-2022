#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate visualization

HG=${BASEDIR}/../genome/hg38.fa

# Plot somatic SVs
for SAMPLE in tumor relapse
do
    bcftools query -f "%CHROM\t%POS\t%CHROM\t%INFO/END\t%ID\t%INFO/SVTYPE\n" ${SAMPLE}.somatic.bcf | egrep "DEL$|INS$|INV$|DUP$" > sv.bed
    bcftools query -f "%CHROM\t%POS\t%INFO/CHR2\t%INFO/POS2\t%ID\t%INFO/SVTYPE\n" ${SAMPLE}.somatic.bcf | grep "BND$" >> sv.bed
    sort -k1,1V -k2,2n sv.bed | cut -f 1-5 | awk '$2>1000 && $4>1000 {print $1"\t"($2-1000)"\t"($2+1000)"\t'${SAMPLE}'_"$5"Left\n"$3"\t"($4-1000)"\t"($4+1000)"\t'${SAMPLE}'_"$5"Right";}' > regions.bed
    rm sv.bed

    # Plot
    wally region -puc -g ${HG} -R regions.bed -s 2 -x 2048 ${BASEDIR}/../alignment/ont/*.bam
    rm regions.bed
done
