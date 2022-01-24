import pandas as pd
from nanoepitools.annotations.enhancers import Enhancers
from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from nanoepitools.annotations.annotations import GFFAnnotationsReader

from mb_analysis.config import module_config
from mb_analysis.reference_cpgs import ReferenceCpGs


"""
Computing some percentages about DMRs for the result section of the paper,
for example the percentage of DMRs that hit promoters, enhancers, etc
"""


def load_samplecomp_hits():
    samplecomp_hits = []
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        samplecomp_hits += [
            {"chrom": line["chromosome"], "start": line["start"], "end": line["end"]}
            for line in pm.read_file(b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, min_diff=min_diff)
        ]
    
    return merge_duplicate_diffmet_hits(samplecomp_hits)


def annotate_dmr_hits(regions, mapped_hits, ref_cpgs):
    for region_key, hits in mapped_hits.items():
        region = regions(region_key)
        for hit in hits:
            start = max(region["start"], hit["start"])
            end = min(region["end"], hit["end"])
            
            hit[f"region_CpGs"] = ref_cpgs.get_CGs(hit["chrom"], start, end, upper=True)
            print(len(hit[f"region_CpGs"]))


def compute_percent_cpgs_in_annotation(regions, all_hits, mapped_hits, ref_cpgs):
    annotate_dmr_hits(regions, mapped_hits, ref_cpgs)
    
    total_cpgs_in_region = len(
        {f"{hit['chrom']}:{cpg}" for hitlist in mapped_hits.values() for hit in hitlist for cpg in hit["region_CpGs"]}
    )
    total_cpgs = len({f"{hit['chrom']}:{cpg}" for hit in all_hits for cpg in hit["CpGs"]})
    
    return total_cpgs_in_region, total_cpgs, total_cpgs_in_region / total_cpgs


def compute_fraction_cpgs_in_cgi(pm, hits):
    
    cgis = pd.read_csv(
        module_config.cgi_regions_file, sep="\t", skiprows=1, names=["chr", "start", "end"], dtype={"chr": str}
    )
    
    cgi_hits = pm.load_regions_hit(cgis, hits=hits, **pm_parameters)
    fraction_cpgs_in_cgi = compute_percent_cpgs_in_annotation(
        lambda x: cgis.loc[x], hits, cgi_hits, ref_cpgs=reference_cpgs
    )
    return fraction_cpgs_in_cgi


def compute_fraction_cpgs_in_promoters(gff, pm, hits):
    promoter_hits = pm.load_promoters_hit(gff, 2000, 500, map_to="transcript", hits=hits)
    
    def promoter_region(transcriptid):
        transcript = gff.get_transcript(transcriptid)
        return {"chrom": transcript.parent.parent.name, "start": transcript.start, "end": transcript.end}
    
    fraction_cpgs_in_promoter = compute_percent_cpgs_in_annotation(
        promoter_region, hits, promoter_hits, ref_cpgs=reference_cpgs
    )
    return fraction_cpgs_in_promoter


def compute_fraction_cpgs_in_enhancers(enhancers, pm, hits):
    enhancers_hit = pm.load_enhancers_hit(enhancers, hits=hits, map_to="index")
    print(len(enhancers_hit))
    
    def enhancer_region(index):
        return enhancers.enhancers_df.loc[index]
    
    fraction_cpgs_in_enhancer = compute_percent_cpgs_in_annotation(
        enhancer_region, hits, enhancers_hit, ref_cpgs=reference_cpgs
    )
    return fraction_cpgs_in_enhancer


if __name__ == "__main__":
    
    min_diff = 0.5
    
    reference_cpgs = ReferenceCpGs()
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    enhancers = Enhancers(enhancers_annotation_file=module_config.enhancer_cerebellum_file)
    enhancers.load()
    enhancers.annotate_nearest_gene(gff)
    enhancers.filter_nearest_gene_none()
    
    pm_parameters = dict(drop_insignificant=True, b_minus_a=True, min_diff=0.5)
    pm_hmm = PycomethOutput(met_comp_file=module_config.pycometh_primary_relapse_file_hmm)
    pm_cgi = PycomethOutput(met_comp_file=module_config.pycometh_primary_relapse_file_cgi)
    
    hits_hmm = list(pm_hmm.read_file(**pm_parameters))
    hits_cgi = list(pm_cgi.read_file(**pm_parameters))
    hits_combined = merge_duplicate_diffmet_hits(hits_hmm + hits_cgi)
    # ugly workaround for a compatibility issue
    for hit in hits_combined:
        hit["chromosome"] = hit["chrom"]
    
    reference_cpgs.add_cpgs_to_dmr_hits(hits_hmm, upper=True)
    reference_cpgs.add_cpgs_to_dmr_hits(hits_cgi, upper=True)
    reference_cpgs.add_cpgs_to_dmr_hits(hits_combined, upper=True)
    
    hmm_cpgs_in_cgi = compute_fraction_cpgs_in_cgi(pm_hmm, hits_hmm)
    cgi_cpgs_in_cgi = compute_fraction_cpgs_in_cgi(pm_cgi, hits_cgi)
    
    cpgs_in_promoter = compute_fraction_cpgs_in_promoters(gff, pm_hmm, hits_combined)
    
    cpgs_in_enhancer = compute_fraction_cpgs_in_enhancers(enhancers, pm_hmm, hits_combined)
