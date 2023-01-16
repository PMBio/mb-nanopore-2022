## Circle-seq analysis

We analyzed Circle-seq short-read paired-end sequencing data using [circle-enrich-finder](https://github.com/henssen-lab/circle-enrich-filter) 1.0.0, to check if TI threads overlap extrachromosomal circular DNA.

Install:
```
# install pgltools
git clone https://github.com/billgreenwald/pgltools.git
# use the v2.2.0 version
export PATH=$PWD/pgltools/sh:$PATH

# create environement
mamba env create -f circleseq.yaml
```


Run the analysis for the PDX primary and relapse:
```
# setup variables
REF={zenodo}/circleseq/annotation/hs37d5_mm10.fa
PRIMARY_SAMPLE=PDXprimary
RELAPSE_SAMPLE=PDXrelapse
PRIMARY_1_FASTQ={zenodo}/circleseq/lane1PrimaryPDX_1_sequence.fq.gz
PRIMARY_2_FASTQ={zenodo}/circleseq/lane1PrimaryPDX_2_sequence.fq.gz
RELAPSE_1_FASTQ={zenodo}/circleseq/lane1RelapsePDX_1_sequence.fq.gz
RELAPSE_2_FASTQ={zenodo}/circleseq/lane1RelapsePDX_2_sequence.fq.gz

# run analysis
bash run.sh
```
