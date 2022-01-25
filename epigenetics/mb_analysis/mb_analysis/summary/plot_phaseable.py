from meth5.meth5 import MetH5File
from mb_analysis.config import module_config
from nanoepitools.plotting.general_plotting import PlotArchiver

import matplotlib.pyplot as plt

if __name__ == "__main__":
    """
    Plotting fraction of reads that were phaseable
    """
    counts = {}
    pa = PlotArchiver("summary", config=module_config)
    samples = ["Germline", "Primary", "Relapse"]
    
    for sample in samples:
        with MetH5File(module_config.meth5_template_file.format(sample=sample), "r") as f:
            counts[sample] = {
                "hp2": (f.h5_fp["reads/read_groups/haplotype"][()] == 2).sum(),
                "hp1": (f.h5_fp["reads/read_groups/haplotype"][()] == 1).sum(),
                "nohp": (f.h5_fp["reads/read_groups/haplotype"][()] == -1).sum(),
            }
    
    plt.figure()
    plt.title("Phaseable reads")
    plt.bar([0, 4, 8], [counts[s]["nohp"] for s in samples], color="gray", label="Unphased")
    plt.bar([1, 5, 9], [counts[s]["hp1"] for s in samples], color="#D87183", label="HP1")
    plt.bar([2, 6, 10], [counts[s]["hp2"] for s in samples], color="#2F618F", label="HP2")
    plt.ylabel("Number of reads")
    plt.xticks([1, 5, 9], labels=samples)
    plt.legend()
    pa.savefig("number_phaseable_reads")
