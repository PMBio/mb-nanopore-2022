from pathlib import Path

import pysam
import pandas as pd
from meth5.meth5 import MetH5File
import tqdm

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
        print(f"Counting {sample}")
        counts[sample] = {}
        tsv_file = Path(f"{module_config.mbdir}2022_reanalysis/whatshap/haplotags_{sample}.tsv")
        """
        
        fastq_dir = Path(f"{module_config.mbdir}2022_reanalysis/mapping/{sample}/")

        for file in tqdm.tqdm(list(fastq_dir.iterdir())):
            if not file.name.endswith(".sorted.filtered.bam"):
                continue
            counts[sample]["total"] += sum(not r.is_supplementary and not r.is_secondary for r in pysam.AlignmentFile(file, "r").fetch())
        """
        phasing = pd.read_csv(tsv_file, sep="\t")
        counts[sample]["hp1"] = (phasing["group"] == "H1").sum()
        counts[sample]["hp2"] = (phasing["group"] == "H2").sum()
        counts[sample]["nohp"] = (phasing["group"] == "none").sum()
        

    def plot_counts(counts, filename):
        pa.figure()
        plt.title("Phaseable reads")
        plt.bar([0, 3, 6], [counts[s]["nohp"] for s in samples], color="#272727", label="Unphased")
        plt.bar([1, 4, 7], [counts[s]["hp1"]+counts[s]["hp2"] for s in samples], color="#90a9b7", label="Phased")
        for x, s in zip([0.5,3.5, 6.5], counts):
            phased = (counts[s]["hp1"]+counts[s]["hp2"])
            ratio = phased / (counts[s]["hp1"]+counts[s]["hp2"]+counts[s]["nohp"])
            plt.text(x, phased * 1.05, f"{ratio*100:.2f}% phased", ha="center")
        plt.ylabel("Number of reads")
        plt.xticks([0.5, 3.5, 6.5], labels=samples)
        plt.ylim(0, 7000000)
        plt.legend()
        pa.savefig(filename)

    plot_counts(counts, "number_phaseable_reads")
    

    # Old values computed from pycoqc output on EBI server
    counts_old = {"Germline": {"hp1": 1720638, "hp2": 1693764, "nohp": 2393480},
                  "Primary":  {"hp1": 3074224, "hp2": 2860592, "nohp": 4217658},
                  "Relapse":  {"hp1": 1821396, "hp2": 1690955, "nohp": 1927170}}

    plot_counts(counts_old, "number_phaseable_reads_old")
