#!/bin/bash

SCRIPT=$(readlink -f "$0")
BASEDIR=$(dirname "$SCRIPT")
export PATH=${BASEDIR}/../conda/bin:${PATH}

source activate assembly

HG=${BASEDIR}/../genome/hg38.fa
THREADS=32

# Assemble, flye
if [ ! -f primary/assembly.fasta ]
then
    flye --out-dir primary --read-error 0.045 --genome-size 2.9g --threads ${THREADS} --nano-hq ../fastq/ont/Primary.fastq.gz
fi

# Assemble, shasta
if [ ! -f ShastaRun/Assembly.fasta ]
then
    curl -O -L https://github.com/chanzuckerberg/shasta/releases/download/0.10.0/shasta-Linux-0.10.0
    chmod ugo+x shasta-Linux-0.10.0
    zcat ../fastq/ont/Primary.fastq.gz > input.fastq
    ./shasta-Linux-0.10.0 --input input.fastq --config Nanopore-May2022
    rm input.fastq
fi

# Align to GRCh38
if [ -f ${HG} ]
then
    for FASTA in primary/assembly.fasta ShastaRun/Assembly.fasta
    do
	ID=`echo ${FASTA} | sed 's/\/.*$//'`
	echo ${ID}
	if [ ! -f ${ID}.hg38.bam ]
	then
	    # Against hg38
	    minimap2 -t ${THREADS} -ax asm5 -L ${HG} ${FASTA} | samtools sort -@ 4 -m 4G > ${ID}.hg38.bam
	    samtools index -@ 4 ${ID}.hg38.bam

	    # Matches
	    alfred bam2match -r ${HG} -o ${ID}.hg38.match.gz ${ID}.hg38.bam
	fi
    done
fi


# Align to assembly
if [ ! -f assembly_lr_mapping.flye.bam ]
then
    minimap2 -ax map-ont -t ${THREADS} primary/assembly.fasta ../fastq/ont/Primary.fastq.gz | samtools sort -@ 4 -m 4G > assembly_lr_mapping.flye.bam
    samtools index -@ 4 assembly_lr_mapping.flye.bam
fi
if [ ! -f assembly_lr_mapping.shasta.bam ]
then
    minimap2 -ax map-ont -t ${THREADS} ShastaRun/Assembly.fasta ../fastq/ont/Primary.fastq.gz | samtools sort -@ 4 -m 4G > assembly_lr_mapping.shasta.bam
    samtools index -@ 4 assembly_lr_mapping.shasta.bam
fi

