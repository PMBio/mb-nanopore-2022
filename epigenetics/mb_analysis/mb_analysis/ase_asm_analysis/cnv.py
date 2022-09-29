import tqdm

from mb_analysis.config import module_config
import pandas as pd
import tqdm
import numpy as np
from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature

def load_gene_cnv(cnv_filename, vcf, gff, replace_chr=True):
    ret_gene_cnv_dict = {}
    with open(cnv_filename, "rt") as fp, tqdm.tqdm() as pbar:
        for line_str in fp:
            if not line_str.startswith("@"):
                break
        cnv_header = line_str.strip().split("\t")
        cur_chrom = None
        
        for line_str in fp:
            pbar.update(1)
            
            line = {k: v for k, v in zip(cnv_header, line_str.strip().split("\t"))}
            line["pos"] = int(line["POSITION"])
            chrom = line["CONTIG"]
            if replace_chr:
                chrom = chrom.replace("chr", "")
            if chrom != cur_chrom:
                cur_chrom = chrom
                gcc_chrom: GFFFeature = gff.chromosomes[chrom]
                seeker = gcc_chrom.get_sorted_in_range_finder()
            
            overlapping_genes = list(seeker.find(start=line["pos"], end=line["pos"] + 1, max_recursion=0))
            if len(overlapping_genes) == 0:
                continue
                
            counts_ref = int(line["REF_COUNT"])
            counts_alt = int(line["ALT_COUNT"])
            try:
                vcf_index = (chrom, line["pos"], line["ALT_NUCLEOTIDE"])
                vcf_row = vcf.vcf_df.loc[vcf_index]
            except:
                continue
            if vcf_row["hp_blood"] == 0:
                continue
            if vcf_row["hp_blood"] == 1:
                counts_hp1 = counts_alt
                counts_hp2 = counts_ref
            elif vcf_row["hp_blood"] == 2:
                counts_hp2 = counts_alt
                counts_hp1 = counts_ref
            
            for overlapping_gene in overlapping_genes:
                gene_id = overlapping_gene.id.split(":")[1]
                counts = ret_gene_cnv_dict.get(gene_id, (0, 0))
                ret_gene_cnv_dict[gene_id] = (counts[0] + counts_hp1, counts[1] + counts_hp2)
    return ret_gene_cnv_dict


def load_cnv_profile(cnv_filename, vcf, chrom, start, end, replace_chr=True, chr_col="CONTIG", pos_col="POSITION", ref_count_col="REF_COUNT", alt_count_col="ALT_COUNT", alt_base_col="ALT_NUCLEOTIDE"):
    ret_y_hp1 = {}
    ret_y_hp2 = {}
    with open(cnv_filename, "rt") as fp, tqdm.tqdm() as pbar:
        for line_str in fp:
            if not line_str.startswith("@"):
                break
        cnv_header = line_str.strip().split("\t")
        
        for line_str in fp:
            pbar.update(1)
            if not line_str.startswith(chrom):
                continue
            line = {k: v for k, v in zip(cnv_header, line_str.strip().split("\t"))}
            line["pos"] = int(line[pos_col])
            chrom = line[chr_col]
            if replace_chr:
                chrom = chrom.replace("chr", "")
                
            if chrom != chrom:
                continue
            if line["pos"] < start:
                continue
            if line["pos"] > end:
                break
                
            counts_ref = int(line[ref_count_col])
            counts_alt = int(line[alt_count_col])
            try:
                vcf_index = (chrom, line["pos"], line[alt_base_col])
                vcf_row = vcf.vcf_df.loc[vcf_index]
            except:
                continue
            if vcf_row["hp_blood"] == 0:
                continue
            if vcf_row["hp_blood"] == 1:
                counts_hp1 = counts_alt
                counts_hp2 = counts_ref
            elif vcf_row["hp_blood"] == 2:
                counts_hp2 = counts_alt
                counts_hp1 = counts_ref
            
            ret_y_hp1[line["pos"]] = ret_y_hp1.get(line["pos"], 0) + counts_hp1
            ret_y_hp2[line["pos"]] = ret_y_hp2.get(line["pos"], 0) + counts_hp2
            
    ret_x = set(ret_y_hp1.keys()).union(set(ret_y_hp2.keys()))
    ret_x = sorted(list(ret_x))
    ret_y_hp1 = [ret_y_hp1.get(pos, 0) for pos in ret_x]
    ret_y_hp2 = [ret_y_hp2.get(pos, 0) for pos in ret_x]
    return np.array(ret_x), np.array(ret_y_hp1), np.array(ret_y_hp2)


def load_tumor_sample_cnv(cnv_filename, gff):
    df = pd.read_csv(cnv_filename, sep="\t")

    cur_chrom = None
    primary = {}
    relapse = {}
    with tqdm.tqdm(total = df.shape[0]) as pbar:
        for _, row in df.iterrows():
            pbar.update(1)
            if row["chr"] != cur_chrom:
                cur_chrom = row["chr"]
                gcc_chrom: GFFFeature = gff.chromosomes[row["chr"]]
                seeker = gcc_chrom.get_sorted_in_range_finder()
        
            overlapping_genes = list(seeker.find(start=row["start"], end=row["end"] + 1, max_recursion=0))
            for gene in overlapping_genes:
                if gene not in primary:
                    primary[gene] = []
                primary[gene].append(row["tumor_CN"])
                
                if gene not in relapse:
                    relapse[gene] = []
                relapse[gene].append(row["relapse_CN"])

    primary_mean = {gene.sanitized_id(): np.mean(c) for gene, c in primary.items()}
    relapse_mean = {gene.sanitized_id(): np.mean(c) for gene, c in relapse.items()}

    ret = pd.DataFrame({"Primary": primary_mean, "Relapse": relapse_mean})
    
    return ret