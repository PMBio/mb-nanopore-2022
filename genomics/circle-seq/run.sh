#!/bin/bash
#SBATCH --job-name=circleseq
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00
#SBATCH --mem=30G

date
hostname

#conda activate circleseq

PRIMARY_SAMPLE=PrimaryPDX
RELAPSE_SAMPLE=RelapsePDX

# mapping
#bash bwa_mem_align.sh $PRIMARY_SAMPLE $REF
#bash bwa_mem_align.sh $RELAPSE_SAMPLE $REF

# circle-enrich-finder
#git clone https://github.com/henssen-lab/circle-enrich-filter.git
# checkout version 1.0.0
cd circle-enrich-filter
git checkout f0fc687c974c9e4308e423f9cbe883344c6b359a
bash run_CircleEnrichFilter.sh -i ${PWD}/../ecDNA/$PRIMARY_SAMPLE/sorted_${PRIMARY_SAMPLE}.dedup.bam -o ${PWD}/../ecDNA/$PRIMARY_SAMPLE/circleEnrich0 -s 0
bash run_CircleEnrichFilter.sh -i ${PWD}/../ecDNA/$RELAPSE_SAMPLE/sorted_${RELAPSE_SAMPLE}.dedup.bam -o ${PWD}/../ecDNA/$RELAPSE_SAMPLE/circleEnrich0 -s 0

#conda deactivate
