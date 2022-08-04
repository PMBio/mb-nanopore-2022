#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly

HG=${BASEDIR}/../genome/hg38.fa
ID=primary
FASTA=primary/assembly.fasta

./calN50.js -L3.1g ${FASTA} > ${ID}.stats.tsv
minimap2 -cxsplice -C5 ${HG} Homo_sapiens.GRCh38.cdna.all.fa > cdna.ref.paf
minimap2 -cxsplice -C5 ${FASTA} Homo_sapiens.GRCh38.cdna.all.fa > ${ID}.asm.paf
paftools.js asmgene -a cdna.ref.paf ${ID}.asm.paf > ${ID}.cdna.tsv
rm ${ID}.asm.paf cdna.ref.paf
