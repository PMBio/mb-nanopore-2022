#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../../conda/bin:${PATH}

source activate assembly

if [ $# -ne 1 ]
then
    echo ""
    echo "Usage: $0 <assembly.fasta>"
    echo ""
    exit -1
fi

ML=11
FILE=${1}
ID=`echo ${FILE} | sed 's/.fa$//'`

mummer -maxmatch -l ${ML} -b -F -c ${FILE} ${FILE}  > ${ID}.matches
python matches.py -m ${ID}.matches  > self.${ID}
Rscript selfalign.R self.${ID}
rm ${ID}.matches self.${ID}

