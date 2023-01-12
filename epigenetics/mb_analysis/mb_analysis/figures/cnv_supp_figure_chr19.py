import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import scipy.stats
from pingouin import partial_corr
from meth5 import MetH5File
import scipy

from nanoepitools.annotations.annotations import GFFAnnotationsReader, GFFFeature
from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.math import fdr_from_pvals
from mb_analysis.ase_asm_analysis.collective_asm import CollectiveAlleleSpecificMethylation
from mb_analysis.ase_asm_analysis.phased_vcf import PhasedVCF
from mb_analysis.ase_asm_analysis.ase import ASE
from mb_analysis.ase_asm_analysis.ase_vs_asm_plotting import plot_ase_vs_asm, plot_ase_vs_asm_nosplit
from mb_analysis.ase_asm_analysis.cnv import load_cnv_profile
from mb_analysis.summary.summarize import get_diffmet_promoter_primary_relapse, get_diffexp
from mb_analysis.summary.outlier_analysis_result import OutlierAnalysisResults
from mb_analysis.config import module_config

if __name__ == "__main__":
    pa = PlotArchiver("supplementary_figures", config=module_config)
    samples = ["Primary", "Relapse"]
    
    vcf_blood = PhasedVCF(module_config.vcf_blood_file)
    vcf_blood.read()

    
    asc =  {sample: load_cnv_profile(module_config.cnv_file_template.format(sample=sample), vcf_blood, "chr19", 0, 60e6, replace_chr=False) for sample in samples}
    

    with pa.open_multipage_pdf("supp_figure_chr19_asc"):
        for sample in samples:
            pa.figure(figsize=(10, 4))
            x,hp1,hp2 = asc[sample]
            y = hp1 / (hp1 + hp2)
            plt.scatter(x, y, s=1)
            plt.ylim(-0.1,1.1)
            pa.savefig()
    