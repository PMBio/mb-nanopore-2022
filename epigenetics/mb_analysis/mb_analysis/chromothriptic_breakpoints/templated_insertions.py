import pandas as pd

from mb_analysis.config import module_config

class TemplatedInsertions:
    """
    Class for loading templated insertion coordinates from a TSV file
    """
    def __init__(self):
        self.templated_insertions_df = None
        pass

    def load(self):
        col_name_type = [
            (0, "r_chr", str),
            (1, "r_start", int),
            (2, "r_end", int),
            (3, "a_contig", str),
            (4, "a_start", int),
            (5, "a_end", int),
        ]

        self.templated_insertions_df = pd.read_csv(
            module_config.templated_insertions_file,
            skiprows=1,
            sep="\t",
            usecols=[c[0] for c in col_name_type],
            names=[c[1] for c in col_name_type],
            dtype={c[1]: c[2] for c in col_name_type},
        )
        self.templated_insertions_df["r_chr"] = self.templated_insertions_df["r_chr"].map(
            lambda x: x.replace("chr", "")
        )

    def get_contigs(self):
        return set(self.templated_insertions_df["a_contig"])

    def __iter__(self):
        return (row for _, row in self.templated_insertions_df.iterrows())
