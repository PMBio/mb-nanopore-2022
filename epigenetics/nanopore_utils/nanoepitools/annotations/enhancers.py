import pandas as pd

def get_nearest_transcript(enhancer, gff):
    nearest_transcript, dist = gff.chromosomes[enhancer["chr"]].get_nearest_feature(
        start=enhancer["start"], end=enhancer["start"], max_recursion=1
    )
    return nearest_transcript, dist


def create_enhancer_label(enhancer, gff):
    if enhancer["chr"] not in gff.chromosomes:
        return "Enhancer"

    nearest_transcript, dist = get_nearest_transcript(enhancer, gff)

    if dist == 0 or dist > 50000:
        return "Enhancer"
    else:
        return f"Enhancer (nearest transcript: {nearest_transcript.name}, {dist} bp away)"


class Enhancers:
    def __init__(self,enhancers_annotation_file):
        self.enhancers_annotation_file = enhancers_annotation_file
    
    def load(self):
        self.enhancers_df = pd.read_csv(
            self.enhancers_annotation_file, sep="\t", usecols=[0, 1, 2], names=["chr", "start", "end"]
        )
        self.enhancers_df["chr"] = self.enhancers_df["chr"].map(lambda x: x.replace("chr", ""))
    
    def annotate_nearest_transcript_label(self, gff):
        self.enhancers_df["label"] = self.enhancers_df.apply(lambda x: create_enhancer_label(x,gff), axis=1)
    
    def annotate_nearest_gene(self, gff, maxdist=1e10):
        def annotate_nearest_gene(enhancer):
            transcript, dist = get_nearest_transcript(enhancer, gff)
            if dist < maxdist:
                return transcript.parent
            else:
                return None
        
        self.enhancers_df["nearest_gene"] = self.enhancers_df.apply(annotate_nearest_gene, axis=1)
    
    def filter_nearest_gene_none(self):
        self.enhancers_df = self.enhancers_df.loc[~pd.isnull(self.enhancers_df["nearest_gene"])]