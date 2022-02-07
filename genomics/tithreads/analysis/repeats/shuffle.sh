#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate variants

if [ $# -ne 4 ]
then
    echo ""
    echo "Usage: $0 [hg19|hg38] <nsim> <source.segments.bed> <features.bed.gz>"
    echo ""
    exit -1
fi


GENOME=${1}
N=${2}
INPUT=${3}
FEATURE=${4}
FEATID=`echo ${FEATURE} | sed 's/^.*\///' | sed 's/.bed.gz$//'`
BREAKID=`echo ${INPUT} | sed 's/^.*\///' | sed 's/.bed$//'`
MOD=`echo "0.1 * ${N}" | bc -l | awk '{print int($1);}'`
TMP=${FEATID}_${BREAKID}

echo "Genome" ${GENOME}
echo "Breakpoints" ${INPUT}
echo "Feature" ${FEATURE}
echo "#Simulations" ${N}

# Shuffle
echo -e "breakpoints\tfeature\tobsexp\tvalue" > shuffle.${TMP}.tsv
zcat ${FEATURE}  | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > feature.${TMP}.bed
cat ${INPUT}  | cut -f 1-3 | sort -k1,1V -k2,2n | uniq > bp.${TMP}.bed
OBS=`bedtools intersect -a bp.${TMP}.bed -b feature.${TMP}.bed -wao | awk '$4!="."' | cut -f 1-3 | sort | uniq | wc -l`
echo "#Observed overlaps" ${OBS}
if [ ${OBS} -gt 1 ]
then
    echo -e "${BREAKID}\t${FEATID}\tobserved\t${OBS}" >> shuffle.${TMP}.tsv
    for i in `seq 1 ${N}`
    do
	if (( $i % ${MOD} == 0 ))
	then
	    echo "Iteration ${i}"
	fi
	VAL=`bedtools shuffle -maxTries 10000 -chrom -noOverlapping -excl ${GENOME}.gaps -i feature.${TMP}.bed -g ${GENOME}.chrom.sizes | bedtools intersect -wao -a bp.${TMP}.bed -b - | awk '$4!="."' | cut -f 1-3 | sort | uniq | wc -l`
	echo -e "${BREAKID}\t${FEATID}\texpected\t${VAL}" >> shuffle.${TMP}.tsv
    done
    # Plot
    Rscript shuffle.R shuffle.${TMP}.tsv
fi
rm feature.${TMP}.bed bp.${TMP}.bed
