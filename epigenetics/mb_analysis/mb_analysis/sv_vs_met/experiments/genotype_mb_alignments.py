import tqdm
from pathlib import Path
from multiprocessing import Queue, Process

import scipy
import numpy as np
import matplotlib.pyplot as plt

from meth5.meth5 import MetH5File
from mb_analysis.config import module_config
from nanoepitools.plotting.general_plotting import PlotArchiver

from mb_analysis.genotyper.genotyper import (
    find_base_in_alignment,
    read_map_ref_alt_other_alignments,
    read_variants,
    BatchedBamFile,
)


def worker_genotype(input_queue: Queue, output_queue: Queue, bam_files, m5_path):
    with BatchedBamFile(bam_files) as batched_bam, MetH5File(m5_path, "r") as mf:
        while True:
            job = input_queue.get()
            if job is None:
                output_queue.put(None)
                break
            chrom, pos, ref, alt = job
            read_map = read_map_ref_alt_other_alignments(batched_bam, chrom, pos, ref, alt)
            num_cat = {k: sum(1 for _, v in read_map.items() if v == k) for k in ["ref", "alt", "other"]}
            if num_cat["ref"] < 5 or num_cat["alt"] < 5 or num_cat["other"] / len(read_map) > 0.25:
                # bad sv
                output_queue.put("SKIP")
                continue
            
            call_group, call_site, call_llr = get_values_around_sv(chrom, pos, read_map, mf)
            diffmet_regions = find_diffmet(call_group, call_site, call_llr)
            
            output_queue.put((chrom, pos, ref, alt, read_map, call_group, call_site, call_llr, diffmet_regions))


def compute_region_difference(call_group, call_site, call_llr, window_start, window_end):
    idx = (call_site >= window_start) & (call_site < window_end)
    window_llr = call_llr[idx]
    window_group = call_group[idx]
    idx_ref = window_group == "ref"
    idx_alt = window_group == "alt"
    if sum(idx_ref) > 10 and sum(idx_alt) > 10:
        bs_ref = (window_llr[idx_ref] > 2).sum() / idx_ref.sum()
        bs_alt = (window_llr[idx_alt] > 2).sum() / idx_alt.sum()
        if abs(bs_ref - bs_alt) > 0.5:
            _, p = scipy.stats.mannwhitneyu(window_llr[idx_ref], window_llr[idx_alt])
            return bs_ref, bs_alt, p
    
    # Insignifcant or low effect size
    return 0, 0, 1


def find_diffmet(call_group, call_site, call_llr, window_size=1000, stride=250, thres=2):
    idx = np.abs(call_llr) > thres
    call_group = call_group[idx]
    call_llr = call_llr[idx]
    call_site = call_site[idx]
    ret = []
    for window_start in range(min(call_site), max(call_site), stride):
        window_end = window_start + window_size
        bs_ref, bs_alt, p = compute_region_difference(call_group, call_site, call_llr, window_start, window_end)
        if p < 0.1:
            if len(ret) > 0:
                if window_start - ret[-1]["end"] < window_size / 2:
                    # continuous region
                    ret[-1]["end"] = window_end
                    continue
            ret.append({"start": window_start, "end": window_end})
    
    for diffmet_region in ret:
        bs_ref, bs_alt, p = compute_region_difference(
            call_group, call_site, call_llr, diffmet_region["start"], diffmet_region["end"]
        )
        diffmet_region["bs_ref"] = bs_ref
        diffmet_region["bs_alt"] = bs_alt
        diffmet_region["bs_diff"] = bs_alt - bs_ref
        diffmet_region["p"] = p
    
    return ret


def get_values_around_sv(chrom, pos, read_map, mf, window_size=500000):
    values_container = mf[chrom].get_values_in_range(pos - window_size, pos + window_size)
    call_group = values_container.get_read_groups(read_group_map=read_map)
    call_site = values_container.get_ranges()[:, 0]
    call_llr = values_container.get_llrs()
    return call_group, call_site, call_llr


def plot_output(call_group, call_site, call_llr, diffmet_regions):
    colors = {"ref": "b", "alt": "magenta"}
    for group in ["ref", "alt"]:
        idx = np.array([g == group for g in call_group])
        plt.scatter(call_site[idx], call_llr[idx], label=str(group), s=1, c=colors[group])
    
    for diffmet_region in diffmet_regions:
        plt.fill_betweenx([-20, 20], diffmet_region["start"], diffmet_region["end"], alpha=0.3, color="r")
    plt.legend()


def worker_output(output_queue: Queue, num_svs, num_input_workers):
    outfile = Path(module_config.sv_vs_met_dir).joinpath("sv_vs_met_diffmet_regions.tsv")
    with pa.open_multipage_pdf("sv_vs_met_llr_window"), open(outfile, "w") as out_f, tqdm.tqdm(total=num_svs) as pbar:
        out_f.write("chrom\tpos\tref\talt\tdiffmet_start\tdiffmet_end\tbs_ref\tbs_alt\tdiffmet_delta_bs\tpval\n")
        while True:
            job = output_queue.get()
            if job is None:
                num_input_workers -= 1
                if num_input_workers == 0:
                    break
                else:
                    continue
            
            pbar.update(1)
            
            if job == "SKIP":
                continue
            
            chrom, pos, ref, alt, read_map, call_group, call_site, call_llr, diffmet_regions = job
            
            if len(diffmet_regions) > 0:
                for diffmet_region in diffmet_regions:
                    out_f.write(
                        "\t".join(
                            [
                                chrom,
                                str(pos),
                                ref,
                                alt,
                                str(diffmet_region["start"]),
                                str(diffmet_region["end"]),
                                f"{diffmet_region['bs_ref']:.2f}",
                                f"{diffmet_region['bs_alt']:.2f}",
                                f"{diffmet_region['bs_diff']:.2f}",
                                f"{diffmet_region['p']:.2f}",
                            ]
                        )
                        + "\n"
                    )
                    out_f.flush()
            pa.figure(figsize=(15, 7))
            plt.title(f"{chrom}:{pos} Ref:{ref} Alt:{alt}")
            plot_output(call_group, call_site, call_llr, diffmet_regions)
            pa.savefig()


if __name__ == "__main__":
    """
    Here we attempt to genotype Nanopore reads based on SOMATIC variants and then write the methylation rates
    so we can investigate whether there are any interesting methylation patterns that correlate with somatic
    variants and not just haplotype:

    Spoiler: Nothing clearly pops up. Most of the time the somatic variants simply correlates with haplotype,
    so using just one sample we can't distinguish whether methylation changes are simply due to paternal/maternal
    alleles or whether they are truly related to the somatic variant.
    """
    
    vcf_file = Path(module_config.vcf_somatic_file)
    svs = read_variants(vcf_file)
    
    sample = "Primary"
    
    bam_dir = Path(module_config.bam_template_dir.format(sample=sample))
    bam_files = [bf for bf in bam_dir.iterdir() if bf.name.endswith(".bam")]
    
    good_svs = 0
    
    sv_list = [(chrom, sv) for chrom in svs for sv in svs[chrom]]
    
    read_group_label_dict = {1: "ref", 2: "alt"}
    read_group_id_dict = {v: k for k, v in read_group_label_dict.items()}
    
    m5_path = module_config.meth5_template_file.format(sample="Primary")
    
    pa = PlotArchiver("sv_vs_met", config=module_config)
    
    workers = 16
    input_queue = Queue(maxsize=workers * 5)
    output_queue = Queue(maxsize=workers * 100)
    
    worker_processes = [
        Process(target=worker_genotype, args=(input_queue, output_queue, bam_files, m5_path)) for _ in range(workers)
    ]
    output_process = Process(target=worker_output, args=(output_queue, len(sv_list), workers))
    
    for p in worker_processes:
        p.start()
    output_process.start()
    
    num_sv_tested = 0
    num_sv_accepted = 0
    for chrom, sv in sv_list:
        num_sv_tested += 1
        # print("Searching for SV: ", chrom, sv)
        pos, ref, alt = sv
        input_queue.put((chrom, pos, ref, alt))
    
    # Deal poison pills
    for _ in worker_processes:
        input_queue.put(None)
    
    output_process.join()
