from typing import Dict

import numpy as np
import tqdm

from meth5.meth5 import MetH5File, MethlyationValuesContainer, compute_betascore

from nanoepitools.annotations.annotations import GFFAnnotationsReader
from mb_analysis.config import module_config
from nanoepitools.math import nangmean, maxabs
from nanoepitools.pycometh_result import PycomethOutput
from nanoepitools.annotations.enhancers import Enhancers


class CollectiveAlleleSpecificMethylation:
    """
    Loads ASM result, and offers functions for mapping to promoters, enhancers by first computing the
    union of all DMRs across all samples and then measuring the differential methylation in each sample for all DMRs.
    """
    
    def __init__(self, gff: GFFAnnotationsReader, samples=module_config.samples):
        self.samples = samples
        self.gff = gff
        self.enhancers = Enhancers(module_config.enhancer_cerebellum_file)
        self.enhancers.load()
        self.enhancers.annotate_nearest_gene(gff, maxdist=3e4)
        self.enhancers.filter_nearest_gene_none()
    
    def __enter__(self):
        self.sample_h5 = {
            sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r") for sample in self.samples
        }
        return self
    
    def __exit__(self, *args):
        for h5 in self.sample_h5.values():
            h5.__exit__(*args)
    
    def is_sufficient_effect(self, bs_diff):
        return np.abs(bs_diff) > 0.5
    
    def get_diffmet_for_region(self, chromosome: str, start: int, end: int, min_sites_per_region=5):
        sample_diff_met = {}
        for sample, h5 in self.sample_h5.items():
            values_container: MethlyationValuesContainer = h5[chromosome].get_values_in_range(start, end)
            hp_betascore = values_container.get_llr_site_readgroup_aggregate(
                group_key="haplotype", aggregation_fun=compute_betascore
            )
            if 1 in hp_betascore and 2 in hp_betascore:
                common_sites = {(x, y) for x, y in hp_betascore[2][1]}.intersection(
                    {(x, y) for x, y in hp_betascore[1][1]}
                )
                mask_common_sites = {hp: [(x, y) in common_sites for x, y in hp_betascore[hp][1]] for hp in {1, 2}}
                site_beta_scores = {hp: hp_betascore[hp][0][mask_common_sites[hp]] for hp in {1, 2}}
                # Remove nans (which occur when llrs are close to 0
                site_beta_scores = {hp: site_beta_scores[hp][~np.isnan(site_beta_scores[hp])] for hp in {1, 2}}
                cur_diff_met = np.nanmean(site_beta_scores[2]) - np.nanmean(site_beta_scores[1])
                if (
                    len(site_beta_scores[1]) >= min_sites_per_region
                    and len(site_beta_scores[1]) >= min_sites_per_region
                    and self.is_sufficient_effect(cur_diff_met)
                ):
                    sample_diff_met[sample] = cur_diff_met
                else:
                    sample_diff_met[sample] = np.nan
            else:
                sample_diff_met[sample] = np.nan
        return sample_diff_met
    
    def read_sample(self, sample, annotation: str, min_diff=0.5):
        pm_file = module_config.pycometh_haplotype_sample_template_file.format(sample=sample)
        print(pm_file)
        pm = PycomethOutput(met_comp_file=pm_file)
        
        if annotation == "promoters":
            gene_hits = pm.load_promoters_hit(
                self.gff, 2000, 500, b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, min_diff=min_diff
            )
        
        elif annotation == "enhancers":
            gene_hits = pm.load_enhancers_hit(
                self.enhancers, b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, min_diff=min_diff
            )
        
        with tqdm.tqdm(total=len(gene_hits)) as pbar:
            for gene, hits in gene_hits.items():
                for hit in hits:
                    sample_diff_met = self.get_diffmet_for_region(hit["chrom"], hit["start"], hit["end"])
                    for s in sample_diff_met:
                        if sample_diff_met[s] != 0 and not np.isnan(sample_diff_met[s]):
                            yield s, gene, sample_diff_met[s]
                pbar.update(1)
    
    def read_all_samples(self, *args, **kwargs) -> Dict[str, Dict[str, float]]:
        sample_gene_diffmet = {s: {} for s in self.samples}
        for sample in self.samples:
            for result_sample, gene, diffmet in self.read_sample(sample, *args, **kwargs):
                if gene not in sample_gene_diffmet[result_sample]:
                    sample_gene_diffmet[result_sample][gene] = []
                sample_gene_diffmet[result_sample][gene].append(diffmet)
        sample_gene_diffmet = {
            s: {gene: maxabs(gene_diffmet) for gene, gene_diffmet in sample_genes.items()}
            for s, sample_genes in sample_gene_diffmet.items()
        }
        return sample_gene_diffmet
