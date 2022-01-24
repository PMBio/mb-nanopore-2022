import numpy as np
from typing import List, Dict
import logging

class AccessibilityEntry:
    def __init__(
        self,
        chrom: str,
        start: int,
        end: int,
        fwd_strand: bool,
        read_name: str,
        llrs: np.ndarray,
    ):
        self.chrom = chrom
        self.read_name = read_name
        self.start = start
        self.end = end
        self.fwd_strand = fwd_strand
        self.llrs = llrs

    def __str__(self):
        code = "_.-°¯"
        curve = "".join([code[int(round(x / 9) + 2)] for x in self.llrs])
        print("%s%s" % ("+" if self.fwd_strand else "-", curve))


def parse_llrs_string(llrs_string):
    llrs = np.array([x for x in llrs_string.encode("ascii")], dtype=float)
    # maps range [33,125] to the range [-20, 20]
    llrs = llrs - 79
    llrs = (np.exp(np.abs(llrs) / 15.0) - 1) * np.sign(llrs)
    return llrs


def parse_line(line):
    line = line.split("\t")
    assert line[1] == "+" or line[1] == "-"
    assert len(line) == 6

    chrom = line[0]
    fwd_strand = line[1] == "+"
    start = int(line[2])
    end = int(line[3])
    read_name = line[4]
    llrs = parse_llrs_string(line[5])

    entry = AccessibilityEntry(
        chrom=chrom,
        start=start,
        end=end,
        read_name=read_name,
        fwd_strand=fwd_strand,
        llrs=llrs,
    )
    return entry


class AccessibilityProfile:
    def __init__(self, entries: List[AccessibilityEntry]):
        self.chrom_dict: Dict[str, List[List[AccessibilityEntry]]] = {}
        for entry in entries:
            if entry.chrom not in self.chrom_dict.keys():
                self.chrom_dict[entry.chrom] = [[], []]  # fwd and rev strand
            strand_index = 0 if entry.fwd_strand else 1
            self.chrom_dict[entry.chrom][strand_index].append(entry)

        for chrom in self.chrom_dict.keys():
            for strand_index in [0, 1]:
                self.chrom_dict[chrom][strand_index] = sorted(
                    self.chrom_dict[chrom][strand_index], key=lambda x: (x.start, x.end)
                )

    @staticmethod
    def read(filename):

        with open(filename, "r") as f:
            # skip header
            f.readline()
            # Parse rest
            entries = [parse_line(x.strip()) for x in f.readlines()]

        ret = AccessibilityProfile(entries)

        return ret

    @staticmethod
    def read_from_filelist(filenames):
        all_entries = []
        for filename in filenames:
            logging.info('Reading %s' % filename)
            with open(filename, "r") as f:
                # skip header
                f.readline()
                # Parse rest
                entries = [parse_line(x.strip()) for x in f.readlines()]
                all_entries = all_entries + entries

        ret = AccessibilityProfile(all_entries)

        return ret
