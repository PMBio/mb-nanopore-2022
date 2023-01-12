#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../../conda/bin:${PATH}

source activate variants

HG=${BASEDIR}/../../genome/hg38.fa

echo "04a5f27a-3b54-4891-83b3-629dba58b227" > reads
/opt/dev/wally/src/wally matches -l -g ${HG} -R reads ${BASEDIR}/../alignment/ont/H021-YD15LC-T1.bam
rm reads
