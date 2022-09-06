#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../../genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

echo -e "id\tstart\tend" > runtime.tsv
for F in ${BASEDIR}/rayas/*/*.log
do
    ID=`echo ${F} | sed 's/.log$//' | sed 's/^.*\///'`
    if [ `cut -f 2 ${BASEDIR}/table.tsv | grep -c -w "${ID}"` -eq 0 ]; then continue; fi
    START=`head -n 1 ${F} | sed 's/\].*//' | sed 's/^\[//'`
    END=`tail -n 1 ${F} | sed 's/\].*//' | sed 's/^\[//'`
    echo -e "${ID}\t${START}\t${END}" >> runtime.tsv
done
Rscript runtime.R
rm runtime.tsv
