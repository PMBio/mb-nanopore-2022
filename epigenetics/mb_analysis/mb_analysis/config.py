class ModuleConfig:
    def __init__(self):
        mbdir = "/omics/groups/OE0540/internal/projects/chromothripsis_medulloblastoma/ont_analysis/"
        self.mbdir = mbdir
        self.pycometh_primary_relapse_file_hmm = (
            f"{mbdir}2022_reanalysis/pycometh/samplecomp/cpg/diffmet/diffmet_600_16_hyp_bs_diff_ihw_yes.tsv"
        )
        self.pycometh_primary_relapse_file_cgi = (
            f"{mbdir}2022_reanalysis/pycometh/samplecomp/cpg/diffmet/CGI_diffmet_hyp_bs_diff_ihw_yes.tsv"
        )
        self.pycometh_primary_relapse_file_promoters = (
            f"{mbdir}2022_reanalysis/pycometh/samplecomp/cpg/diffmet/promoter_based/all_met_comp_adj.tsv"
        )
        
        self.pycometh_primary_relapse_file_ensembl_promoters = (
            f"{mbdir}2022_reanalysis/pycometh/samplecomp/cpg/diffmet/ensembl_promoter_based/all_met_comp_adj.tsv"
        )
        self.pycometh_haplotype_sample_template_file = (
            f"{mbdir}2022_reanalysis//pycometh/asm/cpg/diffmet/{{sample}}_diffmet_600_16_hyp_bs_diff_ihw_yes.tsv"
        )
        self.bam_template_dir = f"{mbdir}2022_reanalysis/mapping/{{sample}}/"
        
        self.meth5_template_file = f"{mbdir}2022_reanalysis//met_merged/{{sample}}_cpg.h5"
        
        self.outlier_analysis_result_file = f"{mbdir}rna/OutlierAnalysis_Medulloblastoma.txt"
        
        self.wasp_ase_file = (
            f"{mbdir}rna/WASP_ase/Run8_promoter_based/CHT_ASM_PromotorBased/cht_results_as.annotated.txt.gz"
        )
        self.wasp_ase_candidate_region_file = f"{mbdir}rna/WASP_ase/promoter_regions_3000_1500.bed"
        self.vcf_blood_file = f"{mbdir}2022_reanalysis/svs/final/blood.phased.vcf"
        self.vcf_somatic_file = f"{mbdir}from_reference/haplotyping/ref_svs/somatic.merged.vcf.gz"
        
        self.grand_summary_file = f"{mbdir}2022_reanalysis/grand_summary.txt"
        self.grand_summary_file_methylation = f"{mbdir}2022_reanalysis/grand_summary_methylation.txt"
        self.rna_wasp_targeted_snps_file = f"{mbdir}rna/haplotyped_snps_near_diffmet_promoter_region.txt"
        self.rna_log_counts_file = f"{mbdir}rna/featureCounts.genes.counts.edgeR.log.txt"
        self.diff_expr_file = f"{mbdir}rna/featureCounts.genes.counts.edgeR.log.txt"
        
        self.fusion_genes_primary_file = f"{mbdir}rna/AS-351452-LR-43985.WGS_guided.fusionGenes.annotated.txt"
        self.fusion_genes_relapse_file = f"{mbdir}rna/Arriba.fusionGenes.AS-351453-LR-43986.complete.annotated.tsv"
        
        self.enhancer_cerebellum_file = f"{mbdir}data/supp/enhanceratlas/Cerebellum.bed"
        self.from_assembly_met_merged_dir = f"{mbdir}2022_reanalysis/from_assembly/met_merged"
        self.templated_insertions_file = f"{mbdir}2022_reanalysis/from_reference/chromthripsis_breakpoints/all_templated_insertions/all_contigs_selected_matches.tsv"
        self.templated_insertions_bed_file_nov2021 = (
            f"{mbdir}from_reference/chromthripsis_breakpoints/templated_insertions_nov2021.bed"
        )
        self.telomere_fusions_insertions_bed_file_dec2021 = (
            f"{mbdir}from_reference/chromthripsis_breakpoints/telomere_fusions.bed"
        )
        self.templated_insertions_meth5_path = f"{mbdir}2022_reanalysis/templated_insertions/met_merged/Primary_cpg.h5"
        self.gff_file = f"{mbdir}data/supp/annotation_chr.gff3"
        self.chaching_dir = f"{mbdir}2022_reanalysis/cache/"
        self.reports_dir = f"{mbdir}2022_reanalysis/reports"
        self.samples = ["Germline", "Primary", "Relapse"]
        self.sample_id_dict_rna = {"Primary": "AS.351452.LR.43985", "Relapse": "AS.351453.LR.43986"}
        self.long_ti_met_merged_dir = f"{mbdir}from_assembly/long_templated_insertions/met_merged"
        self.literature_annotation_dir = f"{mbdir}data/supp/literature"
        self.reference_fasta_file = f"{mbdir}2022_reanalysis/reference/hg38.fa"
        self.reference_softmasked_fasta_file = f"{mbdir}2022_reanalysis/reference/hg38.softmasked.fa"
        self.sv_vs_met_dir = f"{mbdir}2022_reanalysis/sv_vs_met"
        self.asm_vs_ase_dir = f"{mbdir}2022_reanalysis/ase_vs_asm/"
        self.double_minute_m5 = f"{mbdir}2022_reanalysis/double_minute_2parts/met_merged/dm_cpg.h5"
        self.double_minute_twocontigs_parts_path = (
            f"{mbdir}2022_reanalysis/double_minute_2parts/reference/map_to_ref.paf"
        )
        self.methylation450karray_normed_path = (
            f"{mbdir}from_reference/450k_methylation/ILLUMINA450K_Medulo_Beta-ssNOOB.txt.gz"
        )
        self.methylation450karray_metadata_path = f"{mbdir}data/supp/HM450.hg38.manifest.tsv.gz"
        self.methylation850karray_metadata_path = f"{mbdir}data/supp/EPIC.hg38.manifest.tsv.gz"
        
        self.methylation450karray_sample_dict = {"Primary": "X10006823133_R01C01", "Relapse": "X3998523031_R02C01"}
        self.cgi_regions_file = f"{mbdir}2022_reanalysis/pycometh/CGI.bed"
        self.chroms = ["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrM"]
        self.plot_archive_dir = "/home/r933r/snajder/nanoepitools_plots/medulloblastoma_reanalysis"
        self.syncplot_config_file = "/home/r933r/.config/syncplots.yaml"
        
        # Stuff that would need to be updated for rebasecalled analysisy
        self.double_minute_parts_bed_file = f"{mbdir}data/assembly_tobias/double_minute/high_confidence_parts.bed"
        self.cnv_file_template = f"{mbdir}2022_reanalysis/gatk/{{sample}}_asc.tsv"
        self.ct_breakpoints_file = f"{mbdir}from_reference/chromthripsis_breakpoints/complexSVs.primary_tumor.genes.tsv"
        self.cnv_relapse_file = self.cnv_file_template.format(sample="Relapse")
        self.cnv_primary_file = self.cnv_file_template.format(sample="Primary")
        self.goldenpath_cgi_annotation_file = f"{mbdir}data/supp/goldenpath_hg38_cpgIslandExt.txt.gz"
        self.tumor_coverage_file = f"{mbdir}data/cov.10kbp.gz"
    
    def __getitem__(self, key):
        if key in dir(self):
            return getattr(self, key)
    
    def keys(self):
        return self.__dict__.keys()


module_config = ModuleConfig()
