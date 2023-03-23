#!/bin/bash

date
hostname

sample=$1
genome=$2

##get quality of fastqs
mkdir -p ecDNA/${sample}/qc/
fastqc -o ecDNA/${sample}/qc/ ecDNA/${sample}/lane1${sample}_1_sequence.fq.gz ecDNA/${sample}/lane1${sample}_2_sequence.fq.gz

##trim reads
trim_galore -q 20 --phred33 --illumina --paired -j 4 ecDNA/${sample}/lane1${sample}_1_sequence.fq.gz ecDNA/${sample}/lane1${sample}_2_sequence.fq.gz
mv lane1${sample}_1_sequence.fq.gz_trimming_report.txt ecDNA/${sample}/
mv lane1${sample}_1_sequence_val_1.fq.gz ecDNA/${sample}/
mv lane1${sample}_2_sequence.fq.gz_trimming_report.txt ecDNA/${sample}/
mv lane1${sample}_2_sequence_val_2.fq.gz ecDNA/${sample}/

##align with bwa mem
bwa mem -q -t 4 ${genome} ecDNA/${sample}/lane1${sample}_1_sequence_val_1.fq.gz ecDNA/${sample}/lane1${sample}_2_sequence_val_2.fq.gz | samtools sort -@4 -o ecDNA/${sample}/sorted_${sample}.bam -

##index bam file
samtools index ecDNA/${sample}/sorted_${sample}.bam

##mark duplicates
picard MarkDuplicates \
	I=ecDNA/${sample}/sorted_${sample}.bam \
	O=ecDNA/${sample}/sorted_${sample}.dedup.bam \
	M=ecDNA/${sample}/sorted_${sample}.dedup.metrix

##index bam file
samtools index ecDNA/${sample}/sorted_${sample}.dedup.bam

##create coverage overview in bigwig format
bamCoverage --bam ecDNA/${sample}/sorted_${sample}.dedup.bam -o ecDNA/${sample}/sorted_${sample}.dedup.bw

## remove files
rm ecDNA/${sample}/sorted_${sample}.bam*

chmod +x ecDNA/${sample}/*

echo "=========================================================="
echo "Finished on : $(date)"
echo "=========================================================="
