import numpy as np
from mb_analysis.config import module_config
from meth5.meth5 import MetH5File, ChromosomeContainer
from meth5.sparse_matrix import SparseMethylationMatrixContainer


def sample_fun(sample, met_matrix):
    return np.array(
        [f"{sample} (HP{hp.replace('H','')})" if hp != "none" else sample for hp in met_matrix.read_samples]
    )


class MetMatrixLoader:
    def __init__(self):
        self.h5 = {
            sample: MetH5File(module_config.meth5_template_file.format(sample=sample), "r")
            for sample in module_config.samples
        }
    
    def __del__(self):
        for h5file in self.h5.values():
            try:
                h5file.close()
            except:
                pass
    
    def get_merged_matrix(
        self,
        chrom,
        start,
        end,
        with_relapse=True,
        with_germline=True,
        must_overlap_position=None,
        sample_name_fun=sample_fun,
        requires_read_names_for_sample_name_fun=False
    ):
        merged_matrix = None
        for sample in module_config.samples:
            if sample == "Relapse" and not with_relapse:
                continue
            if sample == "Germline" and not with_germline:
                continue
            chrom_container: ChromosomeContainer = self.h5[sample][chrom]
            rg_names = dict(self.h5[sample].h5_fp["reads"]["read_groups"]["haplotype"].attrs)
            met_matrix = chrom_container.get_values_in_range(start, end).to_sparse_methylation_matrix(
                read_groups_key="haplotype", read_read_names=requires_read_names_for_sample_name_fun
            )
            met_matrix.read_samples = np.array([rg_names.get(str(rg), "none") for rg in met_matrix.read_samples])
            
            if must_overlap_position is not None:
                has_over = (
                    np.array(
                        (met_matrix.met_matrix[:, met_matrix.genomic_coord > must_overlap_position] != 0).sum(axis=1)
                    ).flatten()
                    > 0
                )
                
                has_under = (
                    np.array(
                        (met_matrix.met_matrix[:, met_matrix.genomic_coord < must_overlap_position] != 0).sum(axis=1)
                    ).flatten()
                    > 0
                )
                read_idx = has_over & has_under
                met_matrix = met_matrix.get_submatrix_from_read_mask(read_idx)
                if met_matrix.met_matrix.shape[1] == 0:
                    continue
            
            met_matrix.read_samples = sample_name_fun(sample, met_matrix)
            if merged_matrix is None:
                merged_matrix = met_matrix
            else:
                try:
                    merged_matrix = merged_matrix.merge(met_matrix, sample_names="keep")
                except:
                    print("Nothing to merge in")
                    pass
        return merged_matrix
