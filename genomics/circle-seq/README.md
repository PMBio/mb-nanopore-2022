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
REF=../genome/hs37d5_mm10.fa
PRIMARY_SAMPLE=PrimaryPDX
RELAPSE_SAMPLE=RelapsePDX

# run analysis
bash run.sh
```
