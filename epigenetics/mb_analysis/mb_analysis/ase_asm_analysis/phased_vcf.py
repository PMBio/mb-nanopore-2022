import gzip

import pandas as pd


def classify_hp(hp_text):
    if hp_text == "0|1":
        return 2
    elif hp_text == "1|0":
        return 1
    else:
        return 0


class VCFReader:
    def __init__(self, filename):
        self.filename = filename
        self.vcf_fp = None
        self.col_header = None
        self.state = 0
    
    def __enter__(self):
        self.vcf_fp = gzip.open(self.filename, "rt") if self.filename.endswith(".gz") else open(self.filename, "rt")
        return self
    
    def __exit__(self, *args):
        self.vcf_fp.__exit__(*args)
    
    def read_vcf_header(self):
        if self.state != 0:
            raise ValueError("Reading the VCF header must be the first thing to be read")
        
        for line in self.vcf_fp:
            line = line.strip()
            if line[0:2] == "##":
                yield line
            else:
                self.col_header = line
                self.state = 1
                return
    
    def read_columns_header(self):
        if self.state == 0:
            # Skip through VCF header first
            [None for _ in self.read_vcf_header()]
        if self.state > 1:
            raise ValueError("Trying to read the columns header after already having read past it")
        
        self.state = 2
        return self.col_header
    
    def read_content(self):
        print(self.state)
        if self.state < 2:
            # Skip through headers first
            self.read_columns_header()
        
        for line in self.vcf_fp:
            yield line.strip()


class PhasedVCF:
    """
    Reads a VCF with genotype (phasing) information. Only the chromosome, pos, alt, and the genotype are stored.
    """
    
    def __init__(self, vcf_file: str):
        self.vcf_file = vcf_file
    
    def read(self, include_unphased=False):
        with VCFReader(self.vcf_file) as vcf_reader:
            samples = vcf_reader.read_columns_header().split("\t")[9:]
            
            vcf_rows = []
            num_no_gt = 0
            num_unphased = 0
            for line in vcf_reader.read_content():
                if line is None or line == "":
                    break
                
                line = line.split("\t")
                vcf_row = {
                    "chr": line[0],
                    "pos": int(line[1]),
                    "REF": line[3],
                    "ALT": line[4],
                }
                
                try:
                    genotype_column_in_format = line[8].split(":").index("GT")
                except ValueError:
                    num_no_gt += 1
                    # No genotype information for this call
                    continue
                
                vcf_row.update(
                    {
                        f"hp_{sample}": classify_hp(line[9 + i].split(":")[genotype_column_in_format])
                        for i, sample in enumerate(samples)
                    }
                )
                
                if include_unphased or any([vcf_row[f"hp_{sample}"] != 0 for sample in samples]):
                    # Only include it if at least phased in one sample, or if include_unphased flag is set
                    vcf_rows.append(vcf_row)
                else:
                    num_unphased += 1
        
        print(
            f"{len(vcf_rows)} variants loaded, {num_no_gt} variants had no GT tag, "
            f"{num_unphased} variants were unphased"
        )
        
        self.vcf_df = pd.DataFrame(vcf_rows)
        self.vcf_df = self.vcf_df.set_index(["chr", "pos", "ALT"])
        self.vcf_df = self.vcf_df.sort_index()
