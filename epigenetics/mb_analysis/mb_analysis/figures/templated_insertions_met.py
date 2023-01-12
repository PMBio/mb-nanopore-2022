import tqdm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from meth5.meth5 import MetH5File, compute_betascore
from mb_analysis.config import module_config
from nanoepitools.plotting.general_plotting import PlotArchiver


"""
 Very short alignments tend to get a low methylation prediction, which is why we take
 the methylation prediction from two separate analysis. First we use the mapping to reference to
 determine which reads are TI reads and which aren't. Then we compute the methylation rate for non-TI
 reads from these predictions. Secondly, we take the methylation calls for the TIs from a mapping to
 an assembly of templated insertion threads. This way the methylation rates for both groups are
 computed from the longer mappings / better methylation calls.
"""

pa = PlotArchiver("supplementary_figures", config=module_config)

df = pd.read_csv(module_config.templated_insertions_bed_file_nov2021, sep="\t")

f = MetH5File(module_config.meth5_template_file.format(sample="Primary"), "r")

df = df.groupby("clusterid")

with pa.open_multipage_pdf("templated_insertions_methylation_rate"):
    for cluster in df.groups:
        all_ti_readnames = list()
        all_non_ti_bs = []
        for _, row in df.get_group(cluster).iterrows():
            val = f[row["chr"]].get_values_in_range(row["start"], row["end"])
            val_before = f[row["chr"]].get_values_in_range(row["start"] - 2000, row["start"] - 500)
            val_after = f[row["chr"]].get_values_in_range(row["end"] + 500, row["end"] + 2000)
            
            read_names = val.get_read_names()
            reads_middle = val.get_read_ids()
            reads_ti = (
                set(reads_middle).difference(set(val_before.get_read_ids())).difference(set(val_after.get_read_ids()))
            )
            reads_other = set(reads_middle).difference(reads_ti)
            
            idx_ti = np.array([r in reads_ti for r in reads_middle])
            if len(idx_ti) == 0:
                continue
            
            # validation
            read_groups_middle = np.array(val.get_read_groups("haplotype"))
            for hi in set(read_groups_middle):
                idx = read_groups_middle == hi
            
            all_ti_readnames.append(set(read_names[idx_ti]))
            all_non_ti_bs.append(compute_betascore(val.get_llrs()[~idx_ti]))
        
        all_ti_bs_plot = []
        all_non_ti_bs_plot = []
        pa.figure()
        fti = MetH5File(module_config.templated_insertions_meth5_path+".broken", "r")
        for chrom in fti.get_chromosomes():
            print(chrom)
            val = fti[chrom].get_all_values()
            
            if chrom == "tig00000002":
                # Sorry, we had some file corruption there. This is a workaround
                llrs = np.array(list(fti.h5_fp["chromosomes"]["tig00000002"]["llr"][:485]) + [0, 0, 0, 0, 0]+ list(fti.h5_fp["chromosomes"]["tig00000002"]["llr"][490:]))
            else:
                llrs = val.get_llrs()

            print("Getting read names")
            rn_mapping = fti.h5_fp["reads/read_names_mapping"][()]
            read_names = np.array([rn.decode() for rn in rn_mapping[val.get_read_ids()]])
            print("Getting met rates")
            for ti_i, (non_ti_bs, ti_readnames) in tqdm.tqdm(list(enumerate(zip(all_non_ti_bs, all_ti_readnames)))):
                idx_ti = np.array([r in ti_readnames for r in read_names])
                print(idx_ti.sum())
                if idx_ti.sum() > 0:
                    all_ti_bs_plot.append(compute_betascore(llrs[idx_ti]))
                    all_non_ti_bs_plot.append(non_ti_bs)
        
        print("Met diff", np.nanmean(all_ti_bs_plot) - np.nanmean(all_non_ti_bs_plot))
        plt.scatter(all_non_ti_bs_plot, all_ti_bs_plot, c="k")
        plt.title(f"Templated insertion methylation rate, thread {cluster}")
        plt.ylabel("Templated insertion")
        plt.xlabel("Non-Templated insertion")
        plt.plot([0, 1], [0, 1], c="r")
        plt.gca().set_aspect("equal")
        pa.savefig()
