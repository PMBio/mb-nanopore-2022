import pandas as pd

from mb_analysis.config import module_config


class LiftoverAssemblyToRef:
    """
    A function proividing dynamic lift-over of coordinates based on a assembly-to-assembly mapping in PAF format
    """
    
    def __init__(self):
        pass

    def load(self, paf_file, replace_chr=True):
        self.paf_df = pd.read_csv(paf_file,
            sep="\t",
            usecols=[0, 2, 3, 4, 5, 7, 8],
            names=["a_contig", "a_start", "a_end", "dir", "r_chr", "r_start", "r_end"],
            dtype={
                "a_contig": str,
                "a_start": int,
                "a_end": int,
                "dir": str,
                "r_chr": str,
                "r_start": int,
                "r_end": int,
                "qual": int,
            },
        )
        if replace_chr:
            self.paf_df["r_chr"] = self.paf_df["r_chr"].map(lambda x: x.replace("chr", ""))

    def ref_to_assembly(self, chr: str, start: int, end: int):
        r_chr = self.paf_df.groupby("r_chr").get_group(chr)
        r_chr = r_chr.loc[(r_chr["r_start"] < end) & (start < r_chr["r_end"])]
        for _, row in r_chr.iterrows():
            
            if row["dir"] == "+":
                offset_left = max(0, start - row["r_start"])
                yield row["a_contig"], row["a_start"] + offset_left, row["a_start"] + offset_left + end - start
            else:
                offset_right = max(0, row["r_end"] - end)
                yield row["a_contig"], row["a_start"] + offset_right, row["a_start"] + offset_right + end - start
