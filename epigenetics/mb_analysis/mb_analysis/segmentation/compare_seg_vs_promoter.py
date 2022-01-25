from mb_analysis.asm.pm_result_comparer import PMComparer
from nanoepitools.annotations.mapping_utils import MapToPromoter
from mb_analysis.config import module_config

"""
A few statistics comparing differential methylation calling based on methylome segmentation, sequence segmentation (CG),
promoter region (TSS-2000, TSS+500) testing, and ensembl promoter region testing.
"""

if __name__ == "__main__":
    sample_paths = {
        "cgi": module_config.pycometh_primary_relapse_file_cgi,
        "hmm": module_config.pycometh_primary_relapse_file_hmm,
        "promoters": module_config.pycometh_primary_relapse_file_promoters,
        "ensembl_promoters": module_config.pycometh_primary_relapse_file_ensembl_promoters,
    }
    
    pc = PMComparer(sample_paths, "dmr")
    pc.load_hits(min_diff=0.5, drop_insignificant=True, b_minus_a=True)
    pc.load_promoters_hit()
    
    """ Need to use a different sort of mapping for the ensembl genes because we know each one was created from a specific gene"""
    map_to_promoter = MapToPromoter(pc.gff, 2000, 500)
    ensembl_promoters_hit = [
        pc.gff.chromosomes[h["chrom"]]
        .get_nearest_feature(dist_fn=lambda x: map_to_promoter.dist_function(x, h["start"], h["end"]))[0]
        .sanitized_id()
        for h in pc.hits["ensembl_promoters"]
    ]
    ensembl_promoters_hit = {
        gene: [hit for gene_b, hit in zip(ensembl_promoters_hit, pc.hits["ensembl_promoters"]) if gene_b == gene]
        for gene in set(ensembl_promoters_hit)
    }
    pc.promoters_hit["ensembl_promoters"] = ensembl_promoters_hit
    
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
