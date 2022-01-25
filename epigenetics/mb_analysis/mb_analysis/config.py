class ModuleConfig:
    def __init__(self):
        mbdir = "/homes/snajder/data/medulloblastoma/"
        self.pycometh_primary_relapse_file_hmm = (
            "/homes/snajder/data/medulloblastoma/from_reference_incl_supp/pycometh2/readbs_hmm_ihw/all_met_comp_adj.tsv"
        )
        self.pycometh_primary_relapse_file_cgi = (
            "/homes/snajder/data/medulloblastoma/from_reference_incl_supp/pycometh2/readbs_cgi_ihw/all_met_comp_adj.tsv"
        )
        self.pycometh_primary_relapse_file_promoters = (
            "/homes/snajder/data/medulloblastoma/from_reference_incl_supp/pycometh2/readbs_cgi_ihw/all_met_comp_adj"
            ".tsv"
        )
        self.pycometh_primary_relapse_file_ensembl_promoters = "/homes/snajder/data/medulloblastoma/from_reference_incl_supp/pycometh2/readbs_ensembl_promoters/all_met_comp_adj.tsv"
        self.pycometh_haplotype_sample_template_file = f"/homes/snajder/data/medulloblastoma/from_reference_incl_supp/haplotyping/pycometh2/readbs_testall/all_{{sample}}_met_comp_adj.tsv"
        self.bam_template_dir = f"{mbdir}from_reference_incl_supp/mapping/{{sample}}/"
        # self.meth5_template_file = f"{mbdir}from_reference_incl_supp/met_merged/{{sample}}_cpg_rewritten.h5"
        
        self.meth5_template_file = "/homes/snajder/data1/from_reference_incl_supp_met_merged_{sample}_cpg_rewritten.h5"
        self.double_minute_parts_bed_file = (
            f"{mbdir}from_assembly_old/assembly_tobias/double_minute/high_confidence_parts.bed"
        )
        self.outlier_analysis_result_file = f"{mbdir}rna/OutlierAnalysis_Medulloblastoma.txt"
        self.wasp_ase_file = (
            f"{mbdir}rna/WASP_ase/Run8_promoter_based/CHT_ASM_PromotorBased/cht_results_as.annotated.txt.gz"
        )
        self.wasp_ase_candidate_region_file = f"{mbdir}rna/WASP_ase/promoter_regions_3000_1500.bed"
        self.vcf_blood_file = f"{mbdir}from_reference/haplotyping/haplotags_tobias/blood.phased.vcf"
        self.vcf_somatic_file = f"{mbdir}from_reference/haplotyping/ref_svs/somatic.merged.vcf.gz"
        self.cnv_file_template = f"{mbdir}from_reference/haplotyping/cnv/{{sample}}_ase.txt"
        self.cnv_relapse_file = self.cnv_file_template.format(sample="Relapse")
        self.cnv_primary_file = self.cnv_file_template.format(sample="Primary")
        self.grand_summary_file = f"{mbdir}grand_summary.txt"
        self.grand_summary_file_methylation = f"{mbdir}grand_summary_methylation.txt"
        self.rna_wasp_targeted_snps_file = f"{mbdir}rna/haplotyped_snps_near_diffmet_promoter_region.txt"
        self.rna_log_counts_file = f"{mbdir}rna/featureCounts.genes.counts.edgeR.log.txt"
        self.diff_expr_file = f"{mbdir}rna/featureCounts.genes.counts.edgeR.log.txt"
        
        # self.fusion_genes_primary_file = f"{mbdir}rna/Arriba.fusionGenes.AS-351452-LR-43985.complete.annotated.tsv"
        
        self.fusion_genes_primary_file = f"{mbdir}rna/AS-351452-LR-43985.WGS_guided.fusionGenes.annotated.txt"
        self.fusion_genes_relapse_file = f"{mbdir}rna/Arriba.fusionGenes.AS-351453-LR-43986.complete.annotated.tsv"
        
        self.enhancer_cerebellum_file = f"{mbdir}data/supp/enhanceratlas/Cerebellum.bed"
        self.from_assembly_met_merged_dir = f"{mbdir}from_assembly/met_merged"
        self.ct_breakpoints_file = f"{mbdir}from_reference/chromthripsis_breakpoints/complexSVs.primary_tumor.genes.tsv"
        self.templated_insertions_file = (
            f"{mbdir}from_reference/chromthripsis_breakpoints/all_templated_insertions/all_contigs_selected_matches.tsv"
        )
        self.templated_insertions_bed_file_nov2021 = "/homes/snajder/data/medulloblastoma/from_reference/chromthripsis_breakpoints/templated_insertions_nov2021.bed"
        self.telomere_fusions_insertions_bed_file_dec2021 = (
            "/homes/snajder/data/medulloblastoma/from_reference/chromthripsis_breakpoints/telomere_fusions.bed"
        )
        self.templated_insertions_dir = f"{mbdir}from_reference/chromthripsis_breakpoints"
        self.templated_insertions_meth5_path = (
            "/homes/snajder/data/medulloblastoma/from_assembly/long_templated_insertions/met_merged/Primary_cpg.h5"
        )
        self.assembly_to_ref_paf = f"{mbdir}data/assembly_tobias/wtdbg2.genome.ctg_to_ref.paf"
        self.from_reference_met_merged_dir = f"{mbdir}from_reference_incl_supp/met_merged"
        self.from_reference_mapping_dir = f"{mbdir}from_reference_incl_supp/mapping"
        self.gff_file = "/hps/research1/birney/users/adrien/misc/analyses/Medulloblastoma_DNA_promethion/DNA_pipeline_2/results/input/annotation/annotation.gff3"
        self.chaching_dir = f"{mbdir}cache/"
        self.reports_dir = f"{mbdir}from_reference_incl_supp/reports"
        self.segmentation_folder = f"{mbdir}from_reference/haplotyping/segmentation"
        self.met_dir = f"{mbdir}from_reference_incl_supp/met_merged"
        self.samples = ["Germline", "Primary", "Relapse"]
        self.sample_id_dict_rna = {"Primary": "AS.351452.LR.43985", "Relapse": "AS.351453.LR.43986"}
        self.long_ti_met_merged_dir = f"{mbdir}from_assembly/long_templated_insertions" "/met_merged"
        self.literature_annotation_dir = f"{mbdir}data/supp/literature"
        self.reference_fasta_file = "/homes/snajder/data/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa"
        self.sv_vs_met_dir = f"{mbdir}from_reference_incl_supp/sv_vs_met"
        self.asm_vs_ase_dir = f"{mbdir}from_reference_incl_supp/ase_vs_asm/"
        self.double_minute_liftover_m5 = "/homes/snajder/data1/double_minute_met_merged/dm_cpg_liftover.h5"
        self.primary_without_dm_reads_m5 = "/homes/snajder/data1/double_minute_met_merged/primary_without_dm_cpg.h5"
        self.all_cggc_coords_path = (
            "/homes/snajder/data/reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly_cpggpc.pkl.gz"
        )
        self.double_minute_twocontigs_parts_path = (
            "/homes/snajder/data/medulloblastoma/double_minute_2parts/reference/map_to_ref.paf"
        )
        self.double_minute_twocontigs_primary_bam_path = (
            "/homes/snajder/data/medulloblastoma/double_minute_2parts/mapping/Primary.sorted.bam"
        )
        self.double_minute_twocontigs_primary_bam_morefiltered_path = (
            "/homes/snajder/data/medulloblastoma/double_minute_2parts/mapping/Primary.sorted.filtered.bam"
        )
        self.methylation450karray_path = (
            "/homes/snajder/data/medulloblastoma/from_reference/450k_methylation/ILLUMINA450K_Medulo_Beta.txt.gz"
        )
        self.methylation450karray_normed_path = (
            "/homes/snajder/data/medulloblastoma/from_reference/450k_methylation/ILLUMINA450K_Medulo_Beta-ssNOOB.txt.gz"
        )
        self.methylation450karray_metadata_path = (
            "/homes/snajder/data/medulloblastoma/from_reference/450k_methylation/HM450.hg38.manifest.tsv.gz"
        )
        self.methylation850karray_metadata_path = (
            "/homes/snajder/data/reference/resources/infinium_annotation/hg38/EPIC.hg38.manifest.tsv.gz"
        )
        
        self.methylation450karray_sample_dict = {"Primary": "X10006823133_R01C01", "Relapse": "X3998523031_R02C01"}
        self.cgi_regions_file = "/hps/research1/birney/users/adrien/misc/analyses/Medulloblastoma_DNA_promethion/DNA_pipeline_2/results/methylation/pycometh_cgi_finder/CGI.bed"
        self.chroms = [str(i) for i in range(1, 23)] + ["X", "MT"]
        self.plot_archive_dir = "/homes/snajder/data1/plots_medulloblastoma"
    
    def __getitem__(self, key):
        if key in dir(self):
            return getattr(self, key)
    
    def keys(self):
        return self.__dict__.keys()


module_config = ModuleConfig()
