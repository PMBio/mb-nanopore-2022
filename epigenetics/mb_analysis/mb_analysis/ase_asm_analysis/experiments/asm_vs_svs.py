import pysam
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import scipy.stats
from pingouin import partial_corr

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.pycometh_result import PycomethOutput
from nanoepitools.annotations.enhancers import Enhancers
from nanoepitools.math import fdr_from_pvals
from mb_analysis.ase_asm_analysis.collective_asm import CollectiveAlleleSpecificMethylation
from mb_analysis.ase_asm_analysis.phased_vcf import PhasedVCF
from mb_analysis.ase_asm_analysis.ase import ASE
from mb_analysis.ase_asm_analysis.ase_vs_asm_plotting import plot_ase_vs_asm, plot_ase_vs_asm_nosplit
from mb_analysis.ase_asm_analysis.cnv import load_gene_cnv
from mb_analysis.config import module_config

def is_in_range(position_string, chrom, start, end):
    position_string = position_string.split(":")
    if position_string[0] != chrom:
        return False
    return start <= int(position_string[1]) < end
    
    

if __name__ == "__main__":
    somatic_snvs_file = "/homes/snajder/data/medulloblastoma/from_reference/haplotyping/ref_svs/somatic.merged.vcf.gz"
    somatic_svs_file = "/homes/snajder/data/medulloblastoma/from_reference/haplotyping/ref_svs/tumor.somatic.svs.tsv"
    
    pa = PlotArchiver("mb_asm_ase", config=module_config)
    min_abs_diff = 0.25
    search_snv_before = 1000
    search_snv_after = 1000
    
    """ load data """
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()

    enhancers = Enhancers(enhancers_annotation_file=module_config.enhancer_cerebellum_file)
    enhancers.load()
    enhancers.annotate_nearest_gene(gff, maxdist=5000)
    enhancers.filter_nearest_gene_none()

    pm_file = module_config.pycometh_haplotype_sample_template_file.format(sample="Primary")
    pm = PycomethOutput(met_comp_file=pm_file)

    svs = pd.read_csv(somatic_svs_file, sep="\t", names=["bp1", "bp2", "dir1", "dir2", "type", "strandedness"])
    snvs_f = pysam.VariantFile(somatic_snvs_file, "r")
    gene_hits = pm.load_promoters_hit(gff, 2000, 500, b_minus_a=True, drop_insignificant=False, pval_threshold=0.05, min_diff=0.5)

    enhancer_gene_hits = pm.load_enhancers_hit(enhancers, b_minus_a=True, drop_insignificant=False, pval_threshold=0.05, min_diff=0.5)

    snvs_found = set()
    num_hits = 0
    for gene, hits in {**gene_hits, **enhancer_gene_hits}.items():
        for hit in hits:
            chrom = f"chr{hit['chrom']}"
            for variant in snvs_f.fetch(chrom, hit["start"]-search_snv_before, hit["end"]+search_snv_after):
                print("#### FOUND SOMATIC SNV ####")
                print(gene, hit)
                print(variant)
                snvs_found.update({f"{variant.chrom}:{variant.start}"})
            num_hits+=1
            
            for col in ["bp1", "bp2"]:
                for sv in svs[col][svs[col].map(lambda x: is_in_range(x, chrom, hit["start"], hit["end"]))]:
                    print("#### FOUND SOMATIC SV ####")
                    print(gene, hit)
                    print(sv)

import random
def generate_random_hits(num):
    available_chroms = [str(i) for i in range(1,23)]
    for i in range(num):
        random_chrom = random.choice(available_chroms)
        random_start = random.randint(search_snv_before, gff.chromosomes[random_chrom].sorted_children[-1].end)
        random_end = random_start + 500
        yield {"chrom": random_chrom, "start":random_start, "end":random_end}
        
        
random_snvs_found = []
for i in range(100):
    snvs_found = set()
    for hit in generate_random_hits(547):
        chrom = f"chr{hit['chrom']}"
        for variant in snvs_f.fetch(chrom, hit["start"] - search_snv_before, hit["end"] + search_snv_after):
            print("#### FOUND SOMATIC SNV ####")
            print(gene, hit)
            print(variant)
            snvs_found.update({f"{variant.chrom}:{variant.start}"})
        
        for col in ["bp1", "bp2"]:
            for sv in svs[col][svs[col].map(lambda x: is_in_range(x, chrom, hit["start"], hit["end"]))]:
                print("#### FOUND SOMATIC SV ####")
                print(gene, hit)
                print(sv)
    
    random_snvs_found.append(len(snvs_found))