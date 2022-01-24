import numpy as np
import matplotlib.pyplot as plt

from nanoepitools.plotting.general_plotting import PlotArchiver
from nanoepitools.pycometh_result import PycomethOutput, merge_duplicate_diffmet_hits
from mb_analysis.config import module_config


if __name__ == "__main__":
    samplecomp_hits = []
    for pm_file in module_config.pycometh_primary_relapse_file_hmm, module_config.pycometh_primary_relapse_file_cgi:
        pm = PycomethOutput(met_comp_file=pm_file)
        samplecomp_hits += [
            {"chrom": line["chromosome"], "start": line["start"], "end": line["end"], "diff": line["diff"]}
            for line in pm.read_file(b_minus_a=True, drop_insignificant=True, pval_threshold=0.05, min_diff=0.5)
        ]
    
    samplecomp_hits = merge_duplicate_diffmet_hits(samplecomp_hits)
    
    """ Plot overview of DMR segment length """
    pa = PlotArchiver("supplementary_figures", config=module_config)
    
    with pa.open_multipage_pdf("diffmet_segment_length"):
        lengths = np.log10([h["end"] - h["start"] for h in samplecomp_hits])
        pa.figure()
        plt.hist(lengths, bins=200)
        plt.xlabel("Log 10 number of basepairs")
        plt.ylabel("Frequency")
        plt.xlim(0, 5)
        plt.title(f"Primary vs Relapse diffmet segment length")
        pa.savefig()
        
        lengths = np.log10([h["end"] - h["start"] for h in asm_hits["Primary"]])
        pa.figure()
        plt.hist(lengths, bins=200)
        plt.xlabel("Log 10 number of basepairs")
        plt.ylabel("Frequency")
        plt.xlim(0, 5)
        plt.title(f"Primary ASM segment length")
        pa.savefig()
