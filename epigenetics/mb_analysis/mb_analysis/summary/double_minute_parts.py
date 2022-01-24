import pandas as pd
from mb_analysis.config import module_config


class DoubleMinuteParts:
    def __init__(self):
        self.double_minute_df = None
    
    def load(self):
        self.double_minute_df = pd.read_csv(module_config.double_minute_parts_bed_file, sep="\t")
        self.double_minute_df["chr"]  = self.double_minute_df["chr"].map(lambda x: x.replace("chr", ""))
        return self
    
    def __iter__(self):
        return (
            {"chr": row["chr"], "start": row["refstart"], "end": row["refend"]}
            for idx, row in self.double_minute_df.iterrows()
        )
