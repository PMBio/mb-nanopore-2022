#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly

HG=${BASEDIR}/../genome/hg38.fa
THREADS=32

flye --out-dir primary --read-error 0.045 --genome-size 2.9g --threads ${THREADS} --nano-hq ../fastq/ont/Primary.fastq.gz
