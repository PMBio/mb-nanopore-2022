from __future__ import annotations

import numpy as np
import pandas as pd
import re
import scipy.sparse as sp
from pathlib import Path
from typing import List, Union, Dict, Iterable, Tuple, Optional
from meth5.sparse_matrix import SparseMethylationMatrixContainer

import nanoepitools.math as nem


def load_nanopore_metcalls_from_tsv(
    input_folder: Union[str, Path], samples: List[str], mettypes=["cpg", "gpc", "dam"]
) -> Dict[str, Dict[str, pd.DataFrame]]:
    # matches filenames like "footprinting_24_met_cpg.tsv" or
    # "invitro12_2_0_met_dam.tsv"
    filename_regex = re.compile("^(.*)_([0-9]*)_met_([^_]*)\\.tsv$")
    input_folder = Path(input_folder)

    all_sample_met = dict()
    for sample in samples:
        all_sample_met[sample] = dict()
        for mettype in mettypes:
            all_sample_met[sample][mettype] = None

    for f in input_folder.iterdir():
        mitch = filename_regex.match(f.name)
        if mitch is not None:
            sample, batch, mettype = mitch.groups()
            # Only load select sample and mettypes
            if sample not in samples or mettype not in mettypes:
                continue

            met_part = pd.read_csv(f, sep="\t")
            if all_sample_met[sample][mettype] is None:
                all_sample_met[sample][mettype] = met_part
            else:
                all_sample_met[sample][mettype] = all_sample_met[sample][mettype].append(met_part)

    return all_sample_met


def load_merged_nanopore_metcalls(
    input_folder: Union[str, Path], samples: List[str], chroms: List[str], mettypes: List[str] = ["cpg", "gpc", "dam"],
) -> Dict[str, Dict[str, Dict[str, pd.DataFrame]]]:
    """Loads pickled nanopolish methylation calls as pandas dataframes
    and organizes them by sample, chromosome, and methylation type.

    The output is a cascading tree of dictionaries:
         samplename -> (chromosome -> (methylation type -> dataframe))
    :param input_folder: folder containing subfolders for each sample
    :param samples: which samples to include
    :param chroms: which chromosomes to include
    :param mettypes: which methylation types to include (default: cpg, gpc and dam)
    :return: tree of dictionaries with dataframes as leafes
    """
    base_filename = "{chrom}_met_{mettype}.pkl"
    input_folder = Path(input_folder)

    all_sample_met = dict()
    for sample in samples:
        all_sample_met[sample] = dict()

        sample_dir = input_folder.joinpath(sample)

        for chrom in chroms:
            all_sample_met[sample][chrom] = dict()
            for mettype in mettypes:
                all_sample_met[sample][chrom][mettype] = None

                filename = base_filename.format(chrom=chrom, mettype=mettype)
                filepath = sample_dir.joinpath(filename)

                all_sample_met[sample][chrom][mettype] = pd.read_pickle(filepath, compression="gzip")

    return all_sample_met


class RegionFilter:
    def filter(self, allmet):
        return allmet


class IntervalFilter(RegionFilter):
    def __init__(self, start, end):
        self.start = end
        self.end = end

    def filter(self, allmet):
        return allmet.loc[allmet["start"].map(lambda x: self.start <= x < self.end)]


class QuantileFilter(RegionFilter):
    def __init__(self, qfrom, qto):
        self.qfrom = qfrom
        self.qto = qto

    def filter(self, allmet):
        start_locations = list(set(allmet["start"]))
        start = np.quantile(start_locations, self.qfrom)
        end = np.quantile(start_locations, self.qto)
        return allmet.loc[allmet["start"].map(lambda x: start <= x < end)]


class MultiIntervalFilter(RegionFilter):
    def __init__(self, ranges):
        self.ranges = ranges

    def filter(self, allmet):
        passed = allmet["start"].map(lambda x: any([r[0] <= x < r[1] for r in self.ranges]))
        return allmet.loc[passed]


def filter_nanopolish(
    allmet: pd.DataFrame, region_filter: RegionFilter = RegionFilter(), filter_bad_reads: bool = False,
):
    # Apply filter
    allmet = region_filter.filter(allmet)

    # Filtering regions with too high or too low coverage
    # High coverage spikes are likely highly repeated regions (e.g. near
    # telomeres or centromeres)
    if filter_bad_reads:
        print("Filtering bad locations")
        while True:
            # Filter regions that have insanely high coverage
            loc_coverage = allmet[["start", "read_name"]].groupby("start").count()
            loc_coverage_mean = np.mean(loc_coverage["read_name"])
            loc_coverage_std = np.std(loc_coverage["read_name"])
            cov_spike = loc_coverage_mean + loc_coverage_std * 4

            bad_locs = set(loc_coverage.index[(loc_coverage["read_name"] > cov_spike)])
            allmet = allmet[allmet["start"].map(lambda x: x not in bad_locs)]

            # Filter reads that have really low number of sites
            read_coverage = allmet[["start", "read_name"]].groupby("read_name").count()
            bad_reads = set(read_coverage.index[read_coverage["start"] < 10])
            allmet = allmet[allmet["read_name"].map(lambda x: x not in bad_reads)]
            print("New shape:", allmet.shape)
            if len(bad_locs) == 0 and len(bad_reads) == 0:
                return allmet
            if allmet.shape[0] == 0:
                return allmet


def metcall_dataframe_to_llr_matrix(allmet: pd.DataFrame):
    # Detect where reads start and end, and then create sorted list of reads
    allmet = allmet.sort_values(["read_name", "start"])
    read_start = allmet[["read_name", "start"]].groupby("read_name").min().start
    read_start = read_start.sort_values()
    read_names = np.array(read_start.index)

    genomic_coord_start = np.sort(list(set(allmet["start"])))
    genomic_coord_end = np.sort(list(set(allmet["end"])))
    coord_to_index_dict = {genomic_coord_start[i]: i for i in range(len(genomic_coord_start))}
    # Building sparse read vs site methylation matrix
    # As a compromise between memory usage and processing speed,
    # we build dense blocks (fast, but takes memory),
    # then make them sparse (slow, but memory efficient),
    # and then concatenate the sparse blocks
    met_matrix = sp.lil_matrix((len(read_names), len(genomic_coord_start)))
    read_dict = {read_names[i]: i for i in range(len(read_names))}
    cur_rn = ""
    i = 0
    for e in allmet.itertuples():
        if i % 1000000 == 0:
            print("{0:.2f}".format(i / allmet.shape[0] * 100), end=",")
        if i + 1 % 100000000 == 0:
            print()
        i += 1

        if cur_rn != e[5]:
            cur_rn = e[5]
            read_idx = read_dict[cur_rn]
        met_matrix[read_idx, coord_to_index_dict[e[3]]] = e[6]
    met_matrix = sp.csc_matrix(met_matrix)
    return SparseMethylationMatrixContainer(met_matrix, read_names, genomic_coord_start, genomic_coord_end)


def get_only_single_cpg_calls(metcall: pd.DataFrame):
    """Specifically for methylation calls for CpGs. Since Nanopolish
    groups calls that are close to each other, this function filters out
    only single- cpg calls (for things like k-mer uncertainty analysis)

    :param metcall: the dataframe as produced by nanopolish
    :return: Subset of methylation calls which are for single CpG/GpC sites
    (as opposed to grouped calls)
    """
    return metcall.loc[metcall["sequence"].map(lambda x: len(x) == 11)]


def add_sixmer_column(metcall_onecpg: pd.DataFrame):
    """Specifically for methylation calls for CpGs and GpCs. Adds a
    column with the sixmer around the CpG/GpC site. Assumes only single
    CpG/GpC calls.

    :param metcall_onecpg: the dataframe as produced by nanopolish
    """
    assert all(metcall_onecpg["sequence"].map(lambda x: len(x) == 11))
    metcall_onecpg["sixmer"] = metcall_onecpg["sequence"].map(lambda x: x[3:-2])


def compute_kmer_uncertainty(metcall: pd.DataFrame, uncertainty_method="linear"):
    """Computes summary of uncertainties for each k-mer, showing the
    dependence of methylation call uncertainty and sequence context.

    :param metcall: the dataframe as produced by nanopolish.
    :param uncertainty_method: the uncertainty method. See docstring of
    llr_to_uncertainty. Default: 'linear'
    :return: series with kmers as index and uncertainty as values
    """
    metcall_onecpg = get_only_single_cpg_calls(metcall).copy()
    add_sixmer_column(metcall_onecpg)
    metcall_onecpg["uncertainty"] = nem.llr_to_uncertainty(metcall_onecpg["log_lik_ratio"], method=uncertainty_method)

    metcall_onecpg = metcall_onecpg[["sixmer", "uncertainty"]]
    kmer_uncertainty = metcall_onecpg.groupby("sixmer").mean()
    return kmer_uncertainty["uncertainty"]


def count_kmer_incidents(metcall: pd.DataFrame):
    """:param metcall: the dataframe as produced by nanopolis
    :return: series with kmers as index and counts as values."""
    metcall_onecpg = get_only_single_cpg_calls(metcall).copy()
    add_sixmer_column(metcall_onecpg)
    metcall_onecpg = metcall_onecpg[["sixmer", "log_lik_ratio"]]
    kmer_incidents = metcall_onecpg.groupby("sixmer").count()
    return kmer_incidents["log_lik_ratio"]


def compute_kmer_error(metcall: pd.DataFrame, error_method="llr"):
    """:param metcall: the dataframe as produced by nanopolis.

    :param error_method: what method to use to compute the error. Default: 'llr'
    :return: series with kmers as index and errors as values
    """
    metcall_onecpg = get_only_single_cpg_calls(metcall).copy()
    add_sixmer_column(metcall_onecpg)

    if error_method == "llr":
        p = np.abs(metcall_onecpg["log_lik_ratio"])
        p[metcall_onecpg["correct"]] = -p
        metcall_onecpg["error"] = p
    elif error_method == "confusion_rate":
        metcall_onecpg["error"] = ~metcall_onecpg["correct"]
    else:
        raise ValueError("Invalid error method")
    kmer_error = metcall_onecpg.groupby("sixmer").mean()
    return kmer_error["error"]


def compute_read_statistics(
    metcall: pd.DataFrame, compute_bs=False, compute_length=False, llr_threshold=2.5, min_calls=20,
) -> pd.Series:
    """Computes for each read a number of statistics and return a
    dataframe.

    :param metcall: the dataframe as produced by nanopolish
    :param compute_bs: whether to compute beta score
    :param compute_length: whether to compute read length
    :param llr_threshold: exclude calls that are closer than this threshold to
    zero
    :param min_calls: exclude reads with fewer (included) calls
    :return: A pandas series with the read name as index and beta scores
    as values
    """
    to_merge = {}
    if compute_bs:
        metcall_g = metcall[["read_name", "log_lik_ratio"]].copy()
        metcall_g["ismet"] = metcall_g["log_lik_ratio"] > llr_threshold
        metcall_g["isunmet"] = metcall_g["log_lik_ratio"] < -llr_threshold
        metcall_g = metcall_g.groupby("read_name").sum()
        metcall_g["bs"] = metcall_g["ismet"] / (metcall_g["ismet"] + metcall_g["isunmet"])
        to_merge["bs"] = metcall_g.loc[(metcall_g["ismet"] + metcall_g["isunmet"]) > min_calls]["bs"]
    if compute_length:
        metcall_g = metcall[["read_name", "num_motifs"]].copy()
        # metcall_g["len"] =  (metcall_g["end"] - metcall_g["start"] + 1)
        to_merge["length"] = metcall_g.groupby("read_name").sum()["num_motifs"]

    return pd.DataFrame(to_merge)


def compute_read_methylation_betascore(metcall: pd.DataFrame, **kwargs) -> pd.Series:
    """Computes for each read a beta score of methylation.

    :param metcall: the dataframe as produced by nanopolish
    :param llr_threshold: exclude calls that are closer than this threshold to
    zero
    :param min_calls: exclude reads with fewer (included) calls
    :return: A pandas series with the read name as index and beta scores
    as values
    """
    return compute_read_statistics(metcall, compute_bs=True, **kwargs)["bs"]


def aggregate_met_profile(
    met_pairs: Iterable[Tuple[int, float, Optional[int]]], pos_dict: Dict[int, int], window_size: int,
):
    """Computes a methylation profile for a fixed size window in one
    chromosome, given a dictionary that maps genomic coordinates on a
    chromosome to a position in the window. This is a helper function
    for computing an average profile around a type of feature (such as
    transcription start sites)

    :param met_pairs: an iterable of tuples where the first one is the
    genomic position (int) and the second the methylation rate (float). A third optional
    element can be provided for grouped methylation calls. It should contain a list with
    offsets from the start position.
    :param pos_dict: a dictionary that maps genomic positions to
    position in the window
    :param window_size: The size of the window in basepairs
    :returns two numpy arrays of shape (window_size, ). The first one is the totals and
    the second one the count of values for each position in the window.
    """
    if not all([0 <= v < window_size for v in pos_dict.values()]):
        raise ValueError("Trying to map positions outside of window range")

    met_totals = np.zeros(window_size)
    met_counts = np.zeros(window_size)
    for met_row in met_pairs:
        start_pos = met_row[0]
        ratio = met_row[1]
        offsets = [0]
        if len(met_row) == 3:
            # offsets for multiple calls have been provided
            offsets = met_row[2]
        for offset in offsets:
            absolute_pos = start_pos + offset
            if absolute_pos not in pos_dict.keys():
                continue
            relative_pos = pos_dict[absolute_pos]
            met_totals[relative_pos] += ratio
            met_counts[relative_pos] += 1
    return met_totals, met_counts


def aggregate_met_profile_all_chroms(
    met_pairs_per_chrom: Dict[str, Iterable[Tuple[int, float, Optional[int]]]],
    pos_dict_per_chrom: Dict[str, Dict[int, int]],
    window_size: int,
) -> np.ndarray:
    """Computes a methylation profile for a fixed size window over
    multiple chromosomes, given a dictionary per chromosome that maps
    genomic coordinates to a position in the window.

    :param met_pairs_per_chrom: Key is chromosome name, value is iterable of tuples,
    where each tuple consists of the genomic coordinate on that chromosome and the
    methylation rate. A third optional element can be provided for grouped methylation
    calls. It should contain a list with offsets from the start position.
    :param pos_dict_per_chrom: Key is chromosome, value is a dictionary that
    maps the genomic coordinate on that chromosome to a position in the
    :param window_size: The size of the window in basepairs window
    :return: a numpy array of shape (window_size,) containing the average methylation
    rate for each site in the window
    """
    met_totals = np.zeros(window_size)
    met_count = np.zeros(window_size)
    for chrom in met_pairs_per_chrom.keys():
        # Builds a dictionary that maps each genomic position to their relative
        # (to the TSS) position
        pos_dict = pos_dict_per_chrom[chrom]
        met_pairs = met_pairs_per_chrom[chrom]
        met_totals_chrom, met_count_chrom = aggregate_met_profile(met_pairs, pos_dict, window_size)
        met_totals += met_totals_chrom
        met_count += met_count_chrom
    return met_totals, met_count


def compute_average_metrate_profile_all_chrom(*args, **kwargs) -> np.ndarray:
    """Computes a methylation profile for a fixed size window over
    multiple chromosomes, given a dictionary per chromosome that maps
    genomic coordinates to a position in the window.

    :param met_pairs_per_chrom: Key is chromosome name, value is iterable of tuples,
    where each tuple consists of the genomic coordinate on that chromosome and the
    methylation rate. A third optional element can be provided for grouped methylation
    calls. It should contain a list with offsets from the start position.
    :param pos_dict_per_chrom: Key is chromosome, value is a dictionary that
    maps the genomic coordinate on that chromosome to a position in the
    :param window_size: The size of the window in basepairs window
    :return: a numpy array of shape (window_size,) containing the average methylation
    rate for each site in the window
    """
    met_totals, met_count = aggregate_met_profile_all_chroms(*args, **kwargs)
    return met_totals / met_count
