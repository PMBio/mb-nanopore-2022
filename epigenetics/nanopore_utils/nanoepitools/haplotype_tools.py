import pysam
import sys
from pathlib import Path
from typing import List, Union, Dict, Set


class ReadMapping:
    def __init__(
        self,
        read_name: str,
        mapping_contig: str,
        mapping_start: int,
        mapping_end: int,
        haplotype: int,
    ):
        self.read_name = read_name
        self.mapping_contig = mapping_contig
        self.mapping_start = mapping_start
        self.mapping_end = mapping_end
        self.haplotype = haplotype


class PhaseSet:
    def __init__(self, ps_id: int):
        self.ps_id = ps_id
        self.start: int = sys.maxsize
        self.end: int = 0
        self.chrom: str = None
        self.haplotypes: Set[int] = set()
        self.read_mappings: Dict[str, List[ReadMapping]] = {}

    def get_read_mapping(self, read_name: str):
        if read_name in self.read_mappings.keys():
            return self.read_mappings[read_name]
        else:
            return None

    def add_read_mapping(self, mapping: ReadMapping):
        if mapping.read_name not in self.read_mappings.keys():
            self.read_mappings[mapping.read_name] = []
        self.read_mappings[mapping.read_name].append(mapping)
        self.start = min(self.start, mapping.mapping_start)
        self.end = max(self.end, mapping.mapping_end)
        self.haplotypes.add(mapping.haplotype)
        if self.chrom is None:
            self.chrom = mapping.mapping_contig
        else:
            assert self.chrom == mapping.mapping_contig

    def length(self):
        return self.end - self.start

    def get_num_reads_per_haplotype(self):
        num_reads = {h: 0 for h in self.haplotypes}
        for read_mapping_list in self.read_mappings.values():
            for mapping in read_mapping_list:
                num_reads[mapping.haplotype] += 1
        return num_reads

    def get_num_bases_per_haplotype(self):
        num_bases = {h: 0 for h in self.haplotypes}
        for read_mapping_list in self.read_mappings.values():
            for mapping in read_mapping_list:
                num_bases[mapping.haplotype] += (
                    mapping.mapping_end - mapping.mapping_start
                )
        return num_bases


class PhaseSetCollection:
    def __init__(self):
        self.phased: Dict[int, PhaseSet] = {}
        self.unphased: Dict[str, List[ReadMapping]] = {}

    def add_phased_mapping(self, ps_id: int, mapping: ReadMapping):
        if ps_id not in self.phased.keys():
            self.phased[ps_id] = PhaseSet(ps_id)
        self.phased[ps_id].add_read_mapping(mapping)

    def add_unphased_mapping(self, mapping: ReadMapping):
        if mapping.read_name not in self.unphased.keys():
            self.unphased[mapping.read_name] = []
        self.unphased[mapping.read_name].append(mapping)


def extract_read_haplotype_assignment(
    bamfile: Union[str, Path],
    chroms: List[str] = None,
    read_names: List[str] = None,
    return_unphased: bool = False,
) -> PhaseSetCollection:
    phase_sets = PhaseSetCollection()
    with pysam.AlignmentFile(bamfile, "rb") as f:
        segment: pysam.AlignedSegment
        if chroms is None:
            chroms = f.header.references
        for chrom in chroms:
            for segment in f.fetch(contig=chrom):
                # Sanity check: either it has both HP and PS, or neither
                assert (segment.has_tag("HP") and segment.has_tag("PS")) or (
                    not segment.has_tag("HP") and not segment.has_tag("PS")
                )

                if read_names is not None:
                    # If we are interested in specific read names:
                    if segment.query_name not in read_names:
                        continue

                if segment.has_tag("HP"):
                    hp_id = segment.get_tag("HP")
                else:
                    if not return_unphased:
                        continue
                    hp_id = 0

                mapping = ReadMapping(
                    segment.query_name,
                    segment.reference_name,
                    segment.reference_start,
                    segment.reference_end,
                    hp_id,
                )

                if segment.has_tag("HP"):
                    ps_id = "%s_%s" % (chrom, segment.get_tag("PS"))
                    phase_sets.add_phased_mapping(ps_id, mapping)
                else:
                    phase_sets.add_unphased_mapping(mapping)

    return phase_sets
