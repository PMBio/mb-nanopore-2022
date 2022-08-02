SHELL := /bin/bash

# Targets
TARGETS = .conda .mamba .align .variants .phasing .assembly .svintegration .telomere .templins .visualization
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.mamba: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -n base -c conda-forge mamba && touch .mamba

.align: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n align samtools bcftools htslib bwa minimap2 alfred && touch .align

.variants: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n variants samtools bcftools htslib freebayes && touch .variants

.phasing: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n phasing samtools bcftools htslib whatshap cyvcf2 && source activate phasing && pip install cyvcf2 && touch .phasing

.assembly: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n assembly samtools bcftools htslib wtdbg minimap2 alfred flye && touch .assembly

.svintegration: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && mamba create -y -c conda-forge -c bioconda -n sv bcftools samtools delly sniffles minimap2 && source activate sv && python3 -m pip install Truvari && git clone https://github.com/davidebolo1993/VISOR.git && cd VISOR && python setup.py install && cd ../ && pip3 install git+https://github.com/rrwick/Badread.git && touch .svintegration

.telomere: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n telomere samtools htslib minimap2 alfred && touch .telomere

.templins: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n templins samtools htslib minimap2 alfred && touch .templins

.visualization: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n visualization samtools htslib minimap2 wally && touch .visualization

clean:
	rm -f *~

distclean: clean
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/
