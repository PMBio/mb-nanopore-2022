#!/bin/bash

date
hostname

sample=$1
read1=$2
read2=$3
genome=$4

##get quality of fastqs
mkdir -p ${sample}/qc/
fastqc -o ${sample}/qc/ ${read1} ${read2}

##align with bwa mem
bwa mem -q -t 4 ${genome} ${read1} ${read2} | samtools sort -@4 -o sorted_${sample}.bam -

##index bam file
samtools index sorted_${sample}.bam

##create coverage overview in bigwig format
bamCoverage --bam sorted_${sample}.bam -o ${sample}.bw

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
