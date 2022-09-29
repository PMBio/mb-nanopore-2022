import numpy as np
from pyfaidx import Fasta

"""Since in the reanalysis we aligned to an unmasked reference, let's move over the masking from the softmasked
reference. Sanity checks have been performed manually in advance (references are the same length and all bases which
aren't N are the same except for uppercase/lowercase)"""

f_masked = Fasta("/home/r933r/snajder/reference/GRCH38/Homo_sapiens.GRCh38.dna_sm.primary_assembly_chr.fa")
f_unmasked = Fasta("/home/r933r/data/projects/chromothripsis_medulloblastoma/ont_analysis/2022_reanalysis/reference/hg38.fa")

chroms_masked = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY", "chrMT"]
chroms_unmasked = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY", "chrM"]

for chrom_masked, chrom_unmasked in zip(chroms_masked, chroms_unmasked):
    masked = f_masked[chrom_masked][:].seq
    unmasked = f_unmasked[chrom_unmasked][:].seq
    
    new_seq = "".join(u.lower() if m.islower() else u for u, m in zip(unmasked, masked))
    print(f">{chrom_unmasked}")
    print(new_seq)
