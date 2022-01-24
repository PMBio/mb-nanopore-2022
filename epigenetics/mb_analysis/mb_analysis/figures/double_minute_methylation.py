import scipy.stats
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import scipy.stats
from matplotlib.patches import Rectangle, Patch

import matplotlib.pyplot as plt
from nanoepitools.math import fdr_from_pvals
from nanoepitools.plotting.general_plotting import PlotArchiver
from mb_analysis.config import module_config
from meth5.meth5 import MetH5File

"""
Creates double minute methylation figure. At the time of writing this is Figure 4a.
"""

block_config = [
    {"chr": "17", "start": 16030000, "end": 16360000, "id": 1, "color": "#F99B2F"},
    {"chr": "17", "start": 20860000, "end": 22210000, "id": 2, "color": "#E03723"},
    {"chr": "17", "start": 27090000, "end": 27110000, "id": 3, "color": "#9A0262"},
    {"chr": "17", "start": 45140000, "end": 46030000, "id": 4, "color": "#111575"},
    {"chr": "11", "start": 8430000, "end": 8640000, "id": 5, "color": "#126664"},
    {"chr": "11", "start": 15640000, "end": 16180000, "id": 6, "color": "#6FC33D"},
]


def classify_block(chr, start, end):
    for bc in block_config:
        if chr == bc["chr"] and start > bc["start"] and end < bc["end"]:
            return bc["id"]
    return -1


def block_color(id):
    for bc in block_config:
        if bc["id"] == id:
            return bc["color"]


def ordered_args(*args):
    return sorted(list(args))


def is_new_block(new_row, new_row_blockid, old_blocks):
    if len(old_blocks) == 0:
        return True
    old_block = old_blocks[-1]
    if old_block["block"] != new_row_blockid:
        return True
    if (
        np.abs(new_row["ref_end"] - old_block["parts"][-1][0]) > 5000
        and np.abs(new_row["ref_start"] - old_block["parts"][-1][1]) > 5000
    ):
        return True
    return False


def find_reads_spanning_chromothriptic_breakpoints():
    additional_dm_reads = set()
    for contig in contigs:
        row_iter = asm_to_ref.get_group(contig).iterrows()
        _, last_segment = next(row_iter)
        """For each segment of CS11_17"""
        for _, segment in row_iter:
            "Reads in previous segment"
            met_before = mf[last_segment["chrom"]].get_values_in_range(
                last_segment["ref_start"], last_segment["ref_end"]
            )
            read_names_before = set(met_before.get_read_names())
            
            "Reads in next segment"
            met_after = mf[segment["chrom"]].get_values_in_range(segment["ref_start"], segment["ref_end"])
            read_names_after = set(met_after.get_read_names())
            
            " Reads that don't seem to connect the two segments (they don't seem to span the breakpoint)"
            met_cont_before = mf[last_segment["chrom"]].get_values_in_range(
                last_segment["ref_end"], last_segment["ref_end"] + 1000
            )
            read_names_cont_before = set(met_cont_before.get_read_names())
            met_cont_after = mf[segment["chrom"]].get_values_in_range(segment["ref_start"] - 1000, segment["ref_start"])
            read_names_cont_after = set(met_cont_after.get_read_names())
            
            """ Reads that span the breakpoint """
            read_names_shared = read_names_before.intersection(read_names_after)
            read_names_shared = read_names_shared.difference(read_names_cont_before).difference(read_names_cont_after)
            
            additional_dm_reads = additional_dm_reads.union(read_names_shared)
            last_segment = segment
    return additional_dm_reads


def excludenan(x):
    return x[~np.isnan(x)]


def compute_betascore(llrs, llr_threshold=2):
    num_llrs = (np.abs(llrs) > llr_threshold).sum()
    return (llrs > llr_threshold).sum() / num_llrs if num_llrs > 0 else np.nan


def bs_fun(llrs):
    compute_betascore(llrs, llr_threshold=2)


def get_bs_for_region(chrom, start, end):
    dm_reads = set(dm_h5[chrom].get_values_in_range(start, end).get_read_names())
    dm_reads = dm_reads.union(additional_dm_reads)
    
    region_vals = mf[chrom].get_values_in_range(start, end)
    met_matrix = region_vals.to_sparse_methylation_matrix()
    
    normal_reads = set(region_vals.get_read_names())
    normal_reads = normal_reads.difference(dm_reads)
    
    dm_matrix = met_matrix.get_submatrix_from_read_names(dm_reads)
    normal_matrix = met_matrix.get_submatrix_from_read_names(normal_reads)
    
    coords_dm = (dm_matrix.genomic_coord_end + dm_matrix.genomic_coord) / 2
    coords_normal = (normal_matrix.genomic_coord_end + normal_matrix.genomic_coord) / 2
    
    dm_matrix = np.array(dm_matrix.met_matrix.todense())
    normal_matrix = np.array(normal_matrix.met_matrix.todense())
    bs_dm = (dm_matrix > 2).sum(axis=0) / (np.abs(dm_matrix) > 2).sum(axis=0)
    bs_normal = (normal_matrix > 2).sum(axis=0) / (np.abs(normal_matrix) > 2).sum(axis=0)
    
    idx = ~np.isnan(bs_dm)
    bs_dm = bs_dm[idx]
    coords_dm = coords_dm[idx]
    
    idx = ~np.isnan(bs_normal)
    bs_normal = bs_normal[idx]
    coords_normal = coords_normal[idx]
    return coords_dm, bs_dm, coords_normal, bs_normal


def plot_segment_met(contig, bs_interpolation_step=1000):
    cur_plt_x = 0
    normal_color = "#afdde9"
    dm_color = "#de8787"
    
    for _, segment in asm_to_ref.get_group(contig).iterrows():
        chrom = segment["chrom"]
        if chrom not in dm_h5.get_chromosomes() or chrom not in normal_h5.get_chromosomes():
            continue
        
        start = segment["ref_start"]
        end = segment["ref_end"]
        width = end - start
        print("W1: ", width)
        
        coords_dm, bs_dm, coords_normal, bs_normal = get_bs_for_region(chrom, start, end)
        if all([len(x) > 10 for x in (bs_dm, bs_normal)]):
            plt_metrate_x_normal = []
            plt_metrate_y_normal = []
            plt_metrate_x_dm = []
            plt_metrate_y_dm = []
            for coord, bs in zip(coords_dm, bs_dm):
                coord = coord - start + cur_plt_x
                plt_metrate_x_dm.append(coord)
                plt_metrate_y_dm.append(bs)
            
            for coord, bs in zip(coords_normal, bs_normal):
                coord = coord - start + cur_plt_x
                plt_metrate_x_normal.append(coord)
                plt_metrate_y_normal.append(bs)
            
            f = interp1d(plt_metrate_x_normal, plt_metrate_y_normal)
            x = np.arange(min(plt_metrate_x_normal), max(plt_metrate_x_normal), bs_interpolation_step)
            plt.fill_between(x, 0, f(x), color=normal_color, alpha=0.5, linewidth=0)
            f = interp1d(plt_metrate_x_dm, plt_metrate_y_dm)
            x = np.arange(min(plt_metrate_x_dm), max(plt_metrate_x_dm), bs_interpolation_step)
            plt.fill_between(x, 0, f(x), color=dm_color, alpha=0.5, linewidth=0)
        else:
            print("Not enough llrs for: ", classify_block(chrom, start, end))
        
        cur_plt_x += width
    # plt.grid(which="major", axis="y")
    plt.ylabel("Methylation rate")


def plot_segments(contig):
    plt_x_block_patch = []
    cur_plt_x = 0
    for block in blocks[contig]:
        width = sum([p[1] - p[0] for p in block["parts"]])
        color = block_color(block["block"])
        plt_x_block_patch.append(
            Rectangle(
                (cur_plt_x, -0.5),
                width,
                0.40,
                linewidth=1,
                edgecolor="k",
                facecolor=color,
            )
        )
        cur_plt_x += width
    for rect in plt_x_block_patch:
        plt.gca().add_patch(rect)
    plt.ylim(-0.75, 1.1)
    legend_elements = [
        Patch(facecolor=bc["color"], edgecolor="w", linewidth=0, label=f"chr{bc['chr']}:{bc['start']}-{bc['end']}")
        for bc in block_config
    ]
    plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.yticks([0, 0.5, 1], ["0.0", "0.5", "1.0"])
    plt.axhline(y=0, color="k", linewidth=0.5)
    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)
    plt.gca().spines["left"].set_bounds(0, 1)
    plt.gca().spines["bottom"].set_visible(False)
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    plt.tick_params(bottom=True, labelbottom=True)


def plot_diffmet(contig):
    cur_plt_x = 0
    plt_diffmet_x = []
    plt_diffmet_y = []
    plt_diffmet_effect = []
    
    for block in blocks[contig]:
        width = sum([p[1] - p[0] for p in block["parts"]])
        if width < 5000:
            continue
        print("W2: ", width)
        bs_dm = []
        bs_normal = []
        for part in block["parts"]:
            _, bs_dm_part, _, bs_normal_part = get_bs_for_region(block["chrom"], part[0], part[1])
            bs_dm += bs_dm_part.tolist()
            bs_normal += bs_normal_part.tolist()
        bs_dm = np.array(bs_dm)
        bs_normal = np.array(bs_normal)
        stat, p = scipy.stats.mannwhitneyu(bs_dm, bs_normal)
        plt_diffmet_x.append(cur_plt_x + width / 2)
        plt_diffmet_y.append(p)
        effect = np.nanmean(bs_dm) - np.nanmean(bs_normal)
        effect_str = ""
        for diff in np.arange(-1, 0, 0.1):
            if effect < diff:
                # print("Adding - for ", diff)
                effect_str = effect_str + "-"
        for diff in np.arange(0.1, 1, 0.1):
            if effect > diff:
                # print("Adding + for ", diff)
                effect_str = effect_str + "+"
        plt_diffmet_effect.append(effect_str)
        cur_plt_x += width
    plt_diffmet_y = fdr_from_pvals(np.array(plt_diffmet_y))
    for x, p, s in zip(plt_diffmet_x, plt_diffmet_y, plt_diffmet_effect):
        if p < 0.05:
            plt.text(x, -0.7, s, fontsize=16, ha="center")


if __name__ == "__main__":
    asm_to_ref_path = Path(module_config.double_minute_twocontigs_parts_path)
    asm_to_ref = pd.read_csv(
        asm_to_ref_path,
        sep="\t",
        index_col=None,
        usecols=[0, 2, 3, 4, 5, 7, 8, 10, 11],
        names=["ctg", "ctg_start", "ctg_end", "direction", "chrom", "ref_start", "ref_end", "len", "q"],
        dtype={"chrom": str},
    )
    
    contigs = ["ctg1", "ctg2"]
    asm_to_ref = asm_to_ref.sort_values(["ctg", "ctg_start"])
    asm_to_ref = asm_to_ref.groupby("ctg")
    
    pa = PlotArchiver("double_minute", config=module_config)
    
    """
    Merging blocks that are very close (or perhaps just repeats or inversions)
    a little so it's not such a mess of lines
    """
    blocks = {contig: [] for contig in contigs}
    for contig in blocks:
        for i, row in asm_to_ref.get_group(contig).iterrows():
            block = classify_block(row["chrom"], row["ref_start"], row["ref_end"])
            r = [row["ref_start"], row["ref_end"]]
            newblock = is_new_block(row, block, blocks[contig])
            if newblock:
                blocks[contig].append(
                    {
                        "chrom": row["chrom"],
                        "start": row["ctg_start"],
                        "end": row["ctg_end"],
                        "block": block,
                        "parts": [r],
                    }
                )
            else:
                blocks[contig][-1]["end"] = max(row["ctg_end"], blocks[contig][-1]["end"])
                blocks[contig][-1]["start"] = min(row["ctg_start"], blocks[contig][-1]["start"])
                blocks[contig][-1]["parts"].append(r)
    
    mf = MetH5File(module_config.meth5_template_file.format(sample="Primary"), "r")
    additional_dm_reads = find_reads_spanning_chromothriptic_breakpoints()
    dm_reads = set()
    
    dm_h5 = MetH5File(module_config.double_minute_liftover_m5, "r")
    normal_h5 = MetH5File(module_config.primary_without_dm_reads_m5, "r")


with pa.open_multipage_pdf("dm_two_parts_met_rate"):
    for contig in asm_to_ref.groups.keys():
        plt.figure(figsize=(13, 2))
        plt.title(contig)
        plot_segment_met(contig, bs_interpolation_step=500)
        plot_segments(contig)
        plot_diffmet(contig)
        pa.savefig()
