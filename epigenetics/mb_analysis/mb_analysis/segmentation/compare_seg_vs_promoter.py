import tqdm
from mb_analysis.dmr.pm_result_comparer import PMComparer
from nanoepitools.annotations.mapping_utils import MapToPromoter
from mb_analysis.config import module_config
from mb_analysis.reference_cpgs import ReferenceCpGs, format_cpg

"""
A few statistics comparing differential methylation calling based on methylome segmentation, sequence segmentation (CG),
promoter region (TSS-2000, TSS+500) testing, and ensembl promoter region testing.
"""


def map_unique_promoter(hits):
    map_to_promoter = MapToPromoter(pc.gff, 2000, 500)
    promoters_hit = [
        pc.gff.chromosomes[h["chrom"]]
        .get_nearest_feature(dist_fn=lambda x: map_to_promoter.dist_function(x, h["start"], h["end"]))[0]
        .sanitized_id()
        for h in hits
    ]
    promoters_hit = {
        gene: [hit for gene_b, hit in zip(promoters_hit, hits) if gene_b == gene] for gene in set(promoters_hit)
    }
    return promoters_hit


if __name__ == "__main__":
    sample_paths = {
        "cgi": module_config.pycometh_primary_relapse_file_cgi,
        "hmm": module_config.pycometh_primary_relapse_file_hmm,
        "promoters": module_config.pycometh_primary_relapse_file_promoters,
        "ensembl_promoters": module_config.pycometh_primary_relapse_file_ensembl_promoters,
    }
    sample_settings_default = {
        "promoter_before_tss": 2000,
        "promoter_after_tss": 500,
        "b_minus_a": False,
        "drop_insignificant": False,
        "pval_threshold": 0.05,
        "progress": False,
        "min_diff": 0.5,
    }
    sample_settings = {seg: sample_settings_default for seg in sample_paths}
    
    pc = PMComparer(sample_paths, "dmr", sample_settings=sample_settings)
    pc.load_hits(min_diff=0.5, drop_insignificant=False, b_minus_a=True)
    pc.load_promoters_hit()
    
    """ Need to use a different sort of mapping for the ensembl genes because we know each one was created from a specific gene"""
    
    pc.promoters_hit["ensembl_promoters"] = map_unique_promoter(pc.hits["ensembl_promoters"])
    pc.promoters_hit["promoters"] = map_unique_promoter(pc.hits["promoters"])
    
    cgi_plus_hmm_genes = set(pc.promoters_hit["cgi"].keys()).union(set(pc.promoters_hit["hmm"].keys()))
    promoters_genes = set(pc.promoters_hit["promoters"].keys())
    ensembl_promoters_genes = set(pc.promoters_hit["ensembl_promoters"].keys())
    print(f"Testing CGI and HMM combined yields {len(cgi_plus_hmm_genes)} genes")
    print(f"Testing promoter regions yields {len(pc.promoters_hit['promoters'])} genes")
    print(f"Testing ensembl promoter regions yields {len(pc.promoters_hit['ensembl_promoters'])} genes")
    cgi_plus_hmm_minus_promoters = cgi_plus_hmm_genes.difference(promoters_genes)
    print(f"Testing CGI and HMM combined yields {len(cgi_plus_hmm_minus_promoters)} more genes than testing promoters")
    cgi_plus_hmm_minus_ensembl_promoters = cgi_plus_hmm_genes.difference(ensembl_promoters_genes)
    print(
        f"Testing CGI and HMM combined yields {len(cgi_plus_hmm_minus_ensembl_promoters)} genes more than testing ensembl promoters"
    )
    promoters_minus_cgi_plus_hmm = promoters_genes.difference(cgi_plus_hmm_genes)
    print(f"Testing promoters yields {len(promoters_minus_cgi_plus_hmm)} genes not found by CGI+HMM")
    ensembl_promoters_minus_cgi_plus_hmm = ensembl_promoters_genes.difference(cgi_plus_hmm_genes)
    print(f"Testing ensembl promoters yields {len(ensembl_promoters_minus_cgi_plus_hmm)} genes not found by CGI+HMM")

ref_cpgs = ReferenceCpGs()

for seg in pc.hits:
    promoter_hits = list(pc.promoters_hit[seg].values())
    for hits in promoter_hits:
        for hit in hits:
            hit["CpGs"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=True, formatted=True)
    print(
        f"{seg} promoter hits contain {len({cpg for hits in promoter_hits for hit in hits if 'CpGs' in hit for cpg in hit['CpGs']})} CpGs"
    )

for seg in pc.hits:
    for hit in pc.hits[seg]:
        hit["CpGs"] = ref_cpgs.get_CGs(hit["chrom"], hit["start"], hit["end"], upper=True, formatted=True)
    num_cpgs = len({cpg for hit in pc.hits[seg] for cpg in hit['CpGs']})
    print(
        f"{seg} contains a total of {num_cpgs} CpGs"
    )

num_cpgs_hmm_cgi = len({cpg for hit in pc.hits["hmm"] for cpg in hit['CpGs']}.union({cpg for hit in pc.hits["cgi"] for cpg in hit['CpGs']}))
print(
    f"CGI+HMM contain a total of {num_cpgs_hmm_cgi} unique CpGs"
)

promoter_cpgs = set()
for seg in ["cgi", "hmm"]:
    for gene, hits in list(pc.promoters_hit[seg].items()):
        promoter_coords = [(t.start-2000, t.start+500) if t.direction=="+" else (t.end-500, t.start+2000) for t in pc.gff.get_gene(gene).children.values()]
        
        for hit in hits:
            for cpg in hit["CpGs"]:
                cpg_numeric = int(cpg.split(":")[1])
                for start, end in promoter_coords:
                    if start <= cpg_numeric < end:
                        promoter_cpgs.add(cpg)
                        
print(
    f"CGI+HMM promoter hits contain {len(promoter_cpgs)} CpGs"
)

sum_len = sum(hit["end"]-hit["start"] for hit in pc.hits["hmm"]) + sum(hit["end"]-hit["start"] for hit in pc.hits["cgi"])
mean_len = sum_len / (len(pc.hits["hmm"])+len(pc.hits["cgi"]))
print(
    f"CGI+HMM have average length of {mean_len} CpGs"
)