import scipy
import numpy as np
import pandas as pd

from nanoepitools.math import fdr_from_pvals
from mb_analysis.config import module_config


def zscore_to_pval_twosided(z):
    return (1 - scipy.special.ndtr(np.abs(z))) * 2


class OutlierAnalysisResults:
    """
    Support for loading the RNA expression outlier analysis results
    """
    
    def __init__(self):
        self.outlier_df = None
    
    def load(self, filter=True):
        self.outlier_df = pd.read_csv(module_config.outlier_analysis_result_file, sep="\t", index_col="Gene.id")
        # Remove genes for which we have no expression
        if filter:
            self.outlier_df = self.outlier_df.loc[
                self.outlier_df.apply(
                    lambda x: any(
                        x[f"Exp.{sampleid} (log2(tpm))"] > 0 for sampleid in module_config.sample_id_dict_rna.values()
                    ),
                    axis=1,
                )
            ].copy()
        # Remove versioning of gene id
        self.outlier_df.index = self.outlier_df.index.map(lambda x: x.split(".")[0])
        
        for sampleid in module_config.sample_id_dict_rna.values():
            self.outlier_df[f"p_value {sampleid}"] = zscore_to_pval_twosided(self.outlier_df[f"Zscore.{sampleid}"])
            self.outlier_df[f"FDR {sampleid}"] = fdr_from_pvals(self.outlier_df[f"p_value {sampleid}"])
            self.outlier_df[f"direction {sampleid}"] = (
                self.outlier_df[f"Zscore.{sampleid}"].map(lambda x: 0 if np.isnan(x) else np.sign(x)).astype(int)
            )
        return self
    
    def get_genes_outlier(self, sample="Primary", fdr_threshold=0.05):
        sampleid = module_config.sample_id_dict_rna[sample]
        idx = np.abs(self.outlier_df[f"FDR {sampleid}"]) <= fdr_threshold
        return (geneid for geneid in self.outlier_df.loc[idx].index)
    
    def get_outlier_direction(self, sample="Primary", **kwargs):
        sampleid = module_config.sample_id_dict_rna[sample]
        return (
            (gene, self.outlier_df.at[gene, f"direction {sampleid}"])
            for gene in self.get_genes_outlier(sample=sample, **kwargs)
        )
    
    def get_mb_expression(self, sample="Primary"):
        """simply loading the log2 expression for our medulloblastoma samples"""
        sampleid = module_config.sample_id_dict_rna[sample]
        return self.outlier_df[f"Exp.{sampleid} (log2(tpm))"]
