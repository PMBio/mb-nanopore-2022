import numpy as np
from meth5.meth5 import MetH5File

from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from mb_analysis.ase_asm_analysis.collective_asm import CollectiveAlleleSpecificMethylation
from mb_analysis.config import module_config


if __name__ == "__main__":
    
    asm_file = module_config.pycometh_haplotype_sample_template_file.format(sample="Primary")
    pm_asm = PycomethOutput(met_comp_file=asm_file)
    hits = list(pm_asm.read_file(b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, min_diff=0.5))
    hits = merge_duplicate_diffmet_hits(hits)
    
    """ Compute how many DMRs discovered in primary can be confirmed in relapse """
    confirmed = 0
    with MetH5File(module_config.meth5_template_file.format(sample="Relapse"), "r") as f:
        for hit in hits:
            rg_rates = (
                f[hit["chrom"]]
                .get_values_in_range(hit["start"], hit["end"])
                .get_llr_site_readgroup_rate(group_key="haplotype")
            )
            if all([hp in rg_rates for hp in (1, 2)]):
                diff_relapse = np.nanmean(rg_rates[2][0]) - np.nanmean(rg_rates[1][0])
                if abs(diff_relapse) > 0.5 and ((diff_relapse > 0) == (hit["diff"] > 0)):
                    confirmed += 1
    
    print("Primary hits confirmed in relapse: ", confirmed)
