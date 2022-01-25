import gzip
import numpy as np
import tqdm

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.annotations.enhancers import Enhancers

from nanoepitools.annotations.mapping_utils import MapToPromoter


def merge_duplicate_diffmet_hits(oldhits):
    # Merge so we don't count them double
    newhits = []
    a = 0
    b = 0
    for i, hiti in enumerate(oldhits):
        a += 1
        duplicate = -1
        for j, hitj in enumerate(newhits):
            if hiti["start"] < hitj["end"] and hitj["start"] < hiti["end"] and hiti["chrom"] == hitj["chrom"]:
                duplicate = j
                break
        if duplicate >= 0:
            newhits[duplicate] = {
                "chrom": hiti["chrom"],
                "start": min(hiti["start"], newhits[duplicate]["start"]),
                "end": max(hiti["end"], newhits[duplicate]["end"]),
                "diff": np.mean([hiti["diff"], newhits[duplicate]["diff"]]),
            }
        else:
            b += 1
            newhits.append(hiti)
    return newhits


def create_diffmet_entry(pycometh_line):
    diff_met_entry = {}
    diff_met_entry["chrom"] = pycometh_line["chromosome"]
    diff_met_entry["start"] = pycometh_line["start"]
    diff_met_entry["end"] = pycometh_line["end"]
    diff_met_entry["diffmet"] = pycometh_line["diff"]
    return diff_met_entry


class PycomethOutput:
    def __init__(self, met_comp_file):
        self.met_comp_file = met_comp_file
    
    def read_pvals(self):
        with gzip.open(self.met_comp_file, "rt") if self.met_comp_file.endswith(".gz") else open(
            self.met_comp_file, "rt"
        ) as met_comp_fp:
            met_comp_header = met_comp_fp.readline().strip().split("\t")
            for i, line in enumerate(met_comp_fp):
                line = {k: v for k, v in zip(met_comp_header, line.strip().split("\t"))}
                yield float(line["pvalue"])
    
    def read_file(
        self,
        drop_insignificant=True,
        b_minus_a=False,
        retest_fun=None,
        pval_threshold=0.05,
        min_diff=0.25,
        progress=False,
    ):
        N = 0
        if progress:
            with gzip.open(self.met_comp_file, "rt") if self.met_comp_file.endswith(".gz") else open(
                self.met_comp_file, "rt"
            ) as met_comp_fp:
                N = len(met_comp_fp.readlines())
        
        with gzip.open(self.met_comp_file, "rt") if self.met_comp_file.endswith(".gz") else open(
            self.met_comp_file, "rt"
        ) as met_comp_fp:
            met_comp_header = met_comp_fp.readline().strip().split("\t")
            with tqdm.tqdm(disable=(N == 0), total=N) as pbar:
                for i, line in enumerate(met_comp_fp):
                    pbar.update(1)
                    line = {k: v for k, v in zip(met_comp_header, line.strip().split("\t"))}
                    if drop_insignificant and "Significant" not in line["comment"]:
                        continue
                    
                    if retest_fun is not None:
                        recomputed_pval = retest_fun(line["raw_llr_list"])
                        if recomputed_pval > pval_threshold:
                            continue
                    try:
                        line["adj_pvalue"] = float(line["adj_pvalue"])
                        if line["adj_pvalue"] > pval_threshold:
                            continue
                    except:
                        continue
                    
                    line["start"] = int(line["start"])
                    line["end"] = int(line["end"])
                    
                    line["diff"] = eval(line["difference"])
                    if len(line["diff"]) == 0:
                        continue
                    line["diff"] = line["diff"][0]
                    if abs(line["diff"]) < min_diff:
                        continue
                    if not b_minus_a:
                        line["diff"] = -line["diff"]
                    
                    line["chrom"] = line["chromosome"]
                    yield line
    
    def load_promoters_hit(
        self, gff: GFFAnnotationsReader, promoter_before_tss, promoter_after_tss, map_to="gene", hits=None, **kwargs
    ):
        diff_met_table = {}
        map_to_promoter = MapToPromoter(
            gff,
            promoter_before_tss=promoter_before_tss,
            promoter_after_tss=promoter_after_tss,
            input_will_be_sorted=True,
        )
        
        if hits is None:
            hits = list(self.read_file(**kwargs))
        for line in tqdm.tqdm(hits):
            # find promoter
            ids_recorded = set()
            for overlapping_promoter in map_to_promoter(line["chromosome"], line["start"], line["end"]):
                # From transcript to gene
                gene = overlapping_promoter.parent
                # Record every gene just once (if multiple promoters per gene are hit)
                if map_to == "gene":
                    id = gene.sanitized_id()
                elif map_to == "transcript":
                    id = overlapping_promoter.sanitized_id()
                else:
                    raise ValueError("map_to must be gene or transcript")
                
                if id in ids_recorded:
                    continue
                ids_recorded.update({id})
                
                diff_met_list = diff_met_table.get(id, [])
                diff_met_entry = create_diffmet_entry(line)
                diff_met_entry["gene_name"] = gene.name
                diff_met_list.append(diff_met_entry)
                diff_met_table[id] = diff_met_list
        return diff_met_table
    
    def load_gene_bodies_hit(self, gff: GFFAnnotationsReader, **kwargs):
        diff_met_table = {}
        for line in self.read_file(**kwargs):
            gff_chrom: GFFFeature = gff.chromosomes[line["chromosome"]]
            # find promoter
            genes = list(gff_chrom.get_in_range(line["start"], line["end"], max_recursion=0))
            for gene in genes:
                diff_met_list = diff_met_table.get(gene.id, [])
                diff_met_entry = create_diffmet_entry(line)
                diff_met_entry["gene_name"] = gene.name
                diff_met_list.append(diff_met_entry)
                diff_met_table[gene.id] = diff_met_list
        return diff_met_table
    
    def load_enhancers_hit(self, enhancers: Enhancers, map_to="gene", hits=None, **kwargs):
        diff_met_table = {}
        if hits is None:
            hits = self.read_file(**kwargs)
        for line in hits:
            enhancers_hit = enhancers.enhancers_df.loc[
                (enhancers.enhancers_df["chr"] == line["chromosome"])
                & (enhancers.enhancers_df["start"] < line["end"])
                & (line["start"] < enhancers.enhancers_df["end"])
            ]
            
            for index, enhancers_row in enhancers_hit.iterrows():
                
                diff_met_entry = create_diffmet_entry(line)
                if map_to == "gene":
                    gene = enhancers_row["nearest_gene"]
                    diff_met_entry["gene_name"] = gene.name
                    id = gene.sanitized_id()
                elif map_to == "index":
                    id = index
                else:
                    raise ValueError("map_to must be gene or index")
                
                diff_met_list = diff_met_table.get(id, [])
                
                diff_met_list.append(diff_met_entry)
                diff_met_table[id] = diff_met_list
        return diff_met_table
    
    def load_regions_hit(self, regions_df, hits=None, **kwargs):
        diff_met_table = {}
        if hits is None:
            hits = self.read_file(**kwargs)
        for line in hits:
            regions_hit = regions_df.loc[
                (regions_df["chr"] == line["chromosome"])
                & (regions_df["start"] < line["end"])
                & (line["start"] < regions_df["end"])
            ]
            for index, row in regions_hit.iterrows():
                diff_met_list = diff_met_table.get(index, [])
                diff_met_entry = create_diffmet_entry(line)
                diff_met_list.append(diff_met_entry)
                diff_met_table[index] = diff_met_list
        return diff_met_table
