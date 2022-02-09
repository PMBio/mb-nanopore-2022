#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate assembly

for READ in 0645e9e8-3cc3-4f9f-b6f1-382412b524b6
do
     zcat ${BASEDIR}/../../ont/Primary_tumor.fa.gz | grep -w "${READ}" -A 1 > ${READ}.fa
    # Iterate sizes
    for SIZE in 10000 7000 5000 3000 2000 1000 700 500 200
    do
	./subsample_size.sh ${READ} ${SIZE} > read.stats.${SIZE}.txt
    done
done
