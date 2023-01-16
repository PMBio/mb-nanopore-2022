#!/bin/bash

date
hostname

conda activate circleseq

# mapping
bash bwa_mem_align.sh $PRIMARY_SAMPLE $PRIMARY_1_FASTQ $PRIMARY_2_FASTQ $REF
bash bwa_mem_align.sh $PRIMARY_SAMPLE $RELAPSE_1_FASTQ $RELAPSE_2_FASTQ $REF

# circle-enrich-finder
git clone git clone https://github.com/henssen-lab/circle-enrich-filter.git
cd circle-enrich-filter
bash run_CircleEnrichFilter.sh $PRIMARY_SAMPLE/sorted_PDXprimary.dedup.bam $PRIMARY_SAMPLE/circleEnrich0 $REF
bash run_CircleEnrichFilter.sh $RELAPSE_SAMPLE/sorted_PDXprimary.dedup.bam $RELAPSE_SAMPLE/circleEnrich0 $REF

conda deactivate
