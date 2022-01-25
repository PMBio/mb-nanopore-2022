import numpy as np
import scipy.stats

import matplotlib.pyplot as plt
from mb_analysis.config import module_config
from nanoepitools.plotting.general_plotting import PlotArchiver, plot_2d_density
from nanoepitools.annotations.annotations import GFFAnnotationsReader
from mb_analysis.summary.summarize import get_diffmet_promoter_primary_relapse, get_diffexp

"""
Plotting differential expression versus differential methylation for supplementary figure
"""

if __name__ == "__main__":
    pa = PlotArchiver("supplementary_figures", config=module_config)
    
    gff = GFFAnnotationsReader()
    gff.read(module_config.gff_file, only_protein_coding=False)
    gff.build_index()
    
    diffmet = get_diffmet_promoter_primary_relapse(gff, 2000, 500)
    diffexp = get_diffexp()
    
    genes = list(set(diffmet.keys()).intersection(diffexp.keys()))
    diffmet_intersect = np.array([diffmet[k] for k in genes])
    diffexp_intersect = np.array([diffexp[k] for k in genes])
    
    r, p = scipy.stats.pearsonr(diffexp_intersect, diffmet_intersect)
    regression_slope, regression_offset = np.polyfit(diffmet_intersect, diffexp_intersect, 1)
    
    pa.figure()
    plot_2d_density(diffmet_intersect, diffexp_intersect, cmap="jet")
    x = np.arange(-1, 1, 0.01)
    plt.plot(x, regression_slope * x + regression_offset, c="r")
    plt.title("Differential methylation versus differential expression (Relapse - Primary)")
    plt.xlabel("Methylation rate difference")
    plt.ylabel("Log fold change expression")
    pa.savefig("diffmet_diffexp_correlation")
    
    # Highly differentially expressed
    highly_differentially_expressed = diffexp.loc[abs(diffexp) > 5]
    highly_expressed_with_diffmet = list(set(genes).intersection(set(highly_differentially_expressed.index)))
    diffmet_hewdm = np.array([diffmet[k] for k in highly_expressed_with_diffmet])
    diffexp_hewdm = np.array([diffexp[k] for k in highly_expressed_with_diffmet])

    # Highly differentially expressed
    highly_differentially_expressed = diffexp.loc[abs(diffexp) > 2]
    highly_expressed_with_diffmet = list(set(genes).intersection(set(highly_differentially_expressed.index)))
    diffmet_hewdm = np.array([diffmet[k] for k in highly_expressed_with_diffmet])
    diffexp_hewdm = np.array([diffexp[k] for k in highly_expressed_with_diffmet])
    