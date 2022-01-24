from typing import Optional, Dict, List, Tuple, Union

from pysam import AlignmentFile, AlignedSegment


def pysamloader(bam_files: List[str]):
    for f in bam_files:
        yield AlignmentFile(f, "r")


class BatchedBamFile:
    """
    This class loads multiple BAM files using pysam and then provides an interface
    that looks like it's just one bam file.
    """
    
    def __init__(self, bam_files):
        self.bam_files = bam_files
        
        self.bams = list(pysamloader(bam_files))
    
    def __enter__(self):
        for bam in self.bams:
            bam.__enter__()
        return self
    
    def __exit__(self, *args, **kwargs):
        for bam in self.bams:
            bam.__exit__(*args, **kwargs)
    
    def fetch(self, *args, **kwargs):
        for i, bam in enumerate(self.bams):
            for alignment in bam.fetch(*args, **kwargs):
                yield alignment


def read_variants(vcf_file) -> Dict[str, List[Tuple]]:
    variants = {}
    with gzip.open(vcf_file, "rt") as vcf:
        for line in vcf:
            if line.startswith("#"):
                continue
            line = line.split("\t")
            chrom, pos, ref, alt = line[0].replace("chr", ""), int(line[1]), line[3], line[4]
            if chrom not in variants:
                variants[chrom] = []
            variants[chrom].append((pos, ref, alt))
    return variants


def find_base_in_alignment(alignment: AlignedSegment, pos: int) -> Optional[str]:
    idx_q = 0
    idx_r = pos - alignment.reference_start
    # DONT reverse complement it, because minimap already stores reverse complements if the alignment is in reverse
    seq = alignment.query_sequence
    
    for op, l in alignment.cigartuples:
        ref_consumed = op in {0, 2, 3, 7, 8}
        query_consumed = op in {0, 1, 4, 7, 8}
        
        if ref_consumed:
            idx_r -= l
        if query_consumed:
            idx_q += l
        
        if idx_r < 0:
            if query_consumed:
                # base is in query between idx_q-l , idx_q
                base = seq[idx_q + idx_r - 1]
                return base
            else:
                # position has been deleted
                return None


def read_map_ref_alt_other_alignments(
    bam: Union[BatchedBamFile, AlignmentFile], chrom: str, pos: int, ref: str, alt: str
):
    """
    classifies reads as either ref, alt, or other based on a single variant using a bam file
    :param bam: pysam AlignmentFile or BatchedBam file as defined above
    :param chrom: chromsome of the variant to classify based on
    :param pos: chromsome position of the variant to classify base don
    :param ref: reference base/bases.
    :param alt: alternative base/bases.
    :return:
    """
    read_map = {}
    for alignment in bam.fetch(chrom, pos, pos + 1):
        if alignment.query_sequence is None:
            continue
        is_alt = False
        is_ref = False
        is_other = False
        if len(ref) == 1 and len(alt) == 1:
            # SNP
            base = find_base_in_alignment(alignment, pos)
            if base == alt:
                is_alt = True
            elif base == ref:
                is_ref = True
            else:
                is_other = True
        
        if len(ref) > 1 and len(alt) == 1:
            # Deletion
            bases = [find_base_in_alignment(alignment, pos + offset) for offset in range(1, len(alt))]
            if all((base is None for base in bases)):
                is_alt = True
            elif all(bases[i] == ref[i + 1] for i in range(len(bases))):
                is_ref = True
            else:
                is_other = True
        
        if len(ref) == 1 and len(alt) > 1:
            # Insertion
            bases = [find_base_in_alignment(alignment, pos + offset) for offset in range(1, len(alt))]
            if all((base is None for base in bases)):
                is_ref = True
            elif all(bases[i] == alt[i + 1] for i in range(len(bases))):
                is_alt = True
            else:
                is_other = True
        
        if is_ref:
            read_map[alignment.query_name] = "ref"
        if is_alt:
            read_map[alignment.query_name] = "alt"
        if is_other:
            read_map[alignment.query_name] = "other"
    return read_map
