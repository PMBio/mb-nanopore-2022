import tqdm

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature

def load_gene_cnv(cnv_filename, vcf, gff):
    ret_gene_cnv_dict = {}
    with open(cnv_filename, "rt") as fp, tqdm.tqdm() as pbar:
        cnv_header = fp.readline().strip().split("\t")
        cur_chrom = None
        
        for line_str in fp:
            pbar.update(1)
            line = {k: v for k, v in zip(cnv_header, line_str.strip().split("\t"))}
            line["pos"] = int(line["position"])
            chrom = line["contig"].replace("chr", "")
            if chrom != cur_chrom:
                cur_chrom = chrom
                gcc_chrom: GFFFeature = gff.chromosomes[chrom]
                seeker = gcc_chrom.get_sorted_in_range_finder()
            
            overlapping_genes = list(seeker.find(start=line["pos"], end=line["pos"] + 1, max_recursion=0))
            if len(overlapping_genes) == 0:
                continue
                
            counts_ref = int(line["refCount"])
            counts_alt = int(line["altCount"])
            try:
                vcf_index = (chrom, line["pos"], line["altAllele"])
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

