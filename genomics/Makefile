SHELL := /bin/bash

# Targets
TARGETS = .conda .mamba .align .variants .assembly .telomere
PBASE=$(shell pwd)

all:   	$(TARGETS)

.conda:
	wget 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/conda && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.mamba: .conda
	export PATH=${PBASE}/conda/bin:${PATH} && conda install -y -n base -c conda-forge mamba && touch .mamba

.align: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n align samtools bcftools htslib bwa minimap2 alfred && touch .align

.variants: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n variants samtools bcftools htslib delly nanovar sniffles freebayes && touch .variants

.assembly: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n assembly samtools bcftools htslib wtdbg mummer minimap2 alfred && touch .assembly

.telomere: .conda .mamba
	export PATH=${PBASE}/conda/bin:${PATH} && source activate base && mamba create -y -c conda-forge -c bioconda -n telomere samtools htslib minimap2 alfred && touch .telomere

clean:
	rm -f *~

distclean: clean
	rm -rf $(TARGETS) $(TARGETS:=.o) conda/
