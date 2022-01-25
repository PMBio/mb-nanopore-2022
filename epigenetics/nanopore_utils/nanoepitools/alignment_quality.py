from pathlib import Path

import numpy as np
import pysam


def read_alignment(dirname, fileprefix):
    return count_reads_reverse_supplementary(
        [
            f
            for f in Path(dirname).iterdir()
            if f.name.endswith(".bam") and f.name.startswith(fileprefix)
        ]
    )


class AlignmentTypeTotals:
    def __init__(
        self,
        total_single_alignment,
        total_chimeric,
        total_nonchimeric_multi,
        total_inverse_repeats,
        total_only_supp,
    ):
        self.total_single_alignment = total_single_alignment
        self.total_chimeric = total_chimeric
        self.total_nonchimeric_multi = total_nonchimeric_multi
        self.total_inverse_repeats = total_inverse_repeats
        self.total_only_supp = total_only_supp


def classify_alignment(alignments):
    has_fwd_strand = 0
    has_rev_strand = 0
    has_supplementary = 0
    has_primary = 0
    for a in alignments:
        if a[0]:
            has_rev_strand += 1
        else:
            has_fwd_strand += 1
        if a[1]:
            has_supplementary += 1
        else:
            has_primary += 1

    if len(alignments) > 1:
        if has_primary == 0:
            return "only_supp"
        elif has_supplementary == 0:
            return "nonchimeric_multi"
        elif has_fwd_strand > 0 and has_rev_strand > 0:
            return "inverse_repeat"
        else:
            return "chimeric"
    else:
        return "unique"


def categorize_reads(mapping):
    ret = dict()
    for read in mapping.keys():
        alignment_type = classify_alignment(mapping[read])
        if alignment_type not in ret.keys():
            ret[alignment_type] = []
        ret[alignment_type].append(read)
    return ret


def count_types_of_alignments(mapping):
    total_single_alignment = 0
    total_chimeric = 0
    total_nonchimeric_multi = 0
    total_inverse_repeats = 0
    total_only_supp = 0
    for read in mapping.keys():
        alignment_type = classify_alignment(mapping[read])
        if alignment_type == "only_supp":
            total_only_supp += 1
        elif alignment_type == "nonchimeric_multi":
            total_nonchimeric_multi += 1
        elif alignment_type == "inverse_repeat":
            total_inverse_repeats += 1
        elif alignment_type == "chimeric":
            total_chimeric += 1
        elif alignment_type == "unique":
            total_single_alignment += 1
    return AlignmentTypeTotals(
        total_single_alignment,
        total_chimeric,
        total_nonchimeric_multi,
        total_inverse_repeats,
        total_only_supp,
    )


def count_reads_reverse_supplementary(bamfiles, add_qual=False):
    if not isinstance(bamfiles, list):
        bamfiles = [bamfiles]
    read_mapping = dict()
    for bamfile in bamfiles:
        with pysam.AlignmentFile(bamfile, "rb") as f:
            for x in f.fetch():
                if x.query_name not in read_mapping.keys():
                    read_mapping[x.query_name] = list()
                reverse = x.flag & 16 == 16
                supplementary = x.flag & 2048 == 2048
                if add_qual:
                    qual = np.array([c for c in x.qual.encode("ascii")])
                    read_mapping[x.query_name].append((reverse, supplementary, qual))
                else:
                    read_mapping[x.query_name].append((reverse, supplementary))
    return read_mapping
