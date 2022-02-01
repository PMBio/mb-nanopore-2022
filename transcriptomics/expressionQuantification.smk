"""
Snakefile for bulk RNA-seq mapping.

Author: Marc Jan Bonder
Affiliation: DKFZ & EMBL
Study: Medulloblastoma
Date: 01 / 02 / 2022
"""

import glob
import os
from subprocess import run
import pandas as pd
import re
from os.path import join

shell.prefix("set -euo pipefail;") 

def _multi_arg_start(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)

def _multi_arg_end(flag, files):
    flag = " "+flag
    return " ".join(f + flag for f in files)

def _multi_arg_both_ends(flag1, flag2, files):
    flag1 += " "
    flag2 = " "+flag2
    return " ".join(flag1 + f + flag2 for f in files)

def _multi_arg_both_ends_leafcutter(flag1, flag2, files):
    flag1 += " "
    flag2 = " >>"+flag2+";"
    return " ".join(flag1 + f + flag2 for f in files)
#####Here cluster / storage specific flags need to get set ######
## fasta files for reference transcriptome
fasta_unzipped = './Homo_sapiens.GRCh38.dna.primary_assembly.MedulloSample.fa'
fasta_dict = fasta_unzipped.replace('fa', 'dict')
fasta_idx = fasta_unzipped + '.fai'

## define reference file locations
refflat_file = './Homo_sapiens.GRCh38.101.uscdGenePred.txt'
annotation_gtf = './Homo_sapiens.GRCh38.101.gtf'
#apa_bed_file = './annotation/hg19.apadb_v2_processed_bins_UCSC.bed'
#TE_gtf_file = './Homo_sapiens.GRCh37.gene_TE_merged.gtf'
arriba_black_list = './blacklist_hg38_GRCh38_v2.0.0.tsv.gz'
arriba_wgs_info = './tumor_and_relapse_somatic_svs_and_double_minute.tsv' 

## parameters
STAR_GENOME_DIR = './StuttgartSample'
star_genome_files = ['chrLength.txt', 'chrNameLength.txt', 'chrName.txt', 'chrStart.txt', 'exonInfo.tab', 'Genome', 'genomeParameters.txt', 'SA', 'SAindex', 'sjdbInfo.txt', 'sjdbList.fromGTF.out.tab', 'sjdbList.out.tab', 'transcriptInfo.tab']

## define commands
python_cmd = './conda-envs/rnaProcessing/bin/python'
#Star version 2.7.5a
star_cmd = './conda-envs/rnaProcessing/bin/STAR'
#Salmon version 1.3.0
salmon_cmd = './conda-envs/rnaProcessing/bin/salmon'
#Feature counts version (subread 2.0.1)
featurecount_cmd =  './conda-envs/rnaProcessing/bin/featureCounts'
#LeafCutter version:
''
#Picard version: v2.18.7
picard_cmd = './conda-envs/rnaProcessing/bin/picard'
#verifyBam version: v1.1.3
verifyBam_cmd = './conda-envs/rnaProcessing/bin/verifyBamID'
#RnaSeq_QC version: v2.3.5
rnaseqc_cmd = './conda-envs/rnaProcessing/bin/rnaseqc'
#Trim-galore version: v0.6.5
trim_galore_cmd = './conda-envs/rnaProcessing/bin/trim_galore'
##GATK 3.8
gatk_cmd = './conda-envs/rnaProcessing/bin/gatk3'
##Arriba 2
arriba_cmd = './conda-envs/rnaProcessing/bin/arriba'

## parameter objects and samples
CHRS = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22', 'chrM']
RUNS = ['MedulloSample']

excluded_samples = []
dummyVcfStar = './germline.merged_all.het.snp.vcf'
strelkaVariants = './germline.strelka.nochr.multiMerged.sorted.biallelic.snp.vcf.gz'
mergedVariantCalls = './medullo.hg38.wgs_nochr.multiMerged.sorted.biallelic.snp.vcf.gz'

#####Here ends the cluster / storage specific settings ######

SAMPLES_DICT = {}
for run in RUNS:
    tmp = glob.glob("Data/{0}/RAW/*_R1.fastq.gz".format(run))
    tmp =  [os.path.basename(w).replace('_R1.fastq.gz', '') for w in tmp]
    SAMPLES_DICT[run] = [w for w in tmp if w not in excluded_samples]

SAMPLES = glob.glob("Data/*/RAW/*_R1.fastq.gz")
SAMPLES = [os.path.basename(w).replace('_R1.fastq.gz', '') for w in SAMPLES]
SAMPLES = [w for w in SAMPLES if w not in excluded_samples]

#print(SAMPLES)

## targets
star_genome_output = expand('{genome_dir}/{genome_files}', genome_dir=STAR_GENOME_DIR, genome_files=star_genome_files)
trimOut = []
star_bam_output = []
picard_RnaSeqMetrics_output = []
picard_JumpingLibraryMetrics_output = []
verifyBam_output = []
feature_counts_output = []
gatk_ase_output = []
arriba_fusion_gene_output = []

for run in RUNS:
    trimOut.append(expand('Data/{run}/fastq/{sample}_R1_val_1.fq.gz', run=run, sample=SAMPLES_DICT[run]))
    trimOut.append(expand('Data/{run}/fastq/{sample}_R2_val_2.fq.gz', run=run, sample=SAMPLES_DICT[run]))
    star_bam_output.append(expand('Data/{run}/star/{sample}/{sample}.Aligned.out.bam', run=run, sample=SAMPLES_DICT[run]))
    gatk_ase_output.append(expand('Data/{run}/ase/{sample}.ase.all.tsv', run=run, sample=SAMPLES_DICT[run]))
    gatk_ase_output.append(expand('Data/{run}/ase/{sample}.ase.strelka.tsv', run=run, sample=SAMPLES_DICT[run]))
    arriba_fusion_gene_output.append(expand('Data/{run}/arriba/{sample}.fusionGenes.tsv', run=run, sample=SAMPLES_DICT[run]))
    arriba_fusion_gene_output.append(expand('Data/{run}/arriba/{sample}.WGS_guided.fusionGenes.tsv', run=run, sample=SAMPLES_DICT[run]))
    arriba_fusion_gene_output.append(expand('Data/{run}/arriba/{sample}.WGS_guided.noBlackList.fusionGenes.tsv', run=run, sample=SAMPLES_DICT[run]))
    picard_RnaSeqMetrics_output.append(expand('Data/{run}/star/{sample}/{sample}.Aligned.rna_metrics', run=run, sample=SAMPLES_DICT[run]))
    picard_JumpingLibraryMetrics_output.append(expand('Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.jump_metrics', run=run, sample=SAMPLES_DICT[run]))
    verifyBam_output.append(expand('Data/{run}/verifyBamId/{sample}.bestSM', run=run, sample=SAMPLES_DICT[run]))
    feature_counts_output.append(expand('Data/{run}/initialQuant_featureCounts/{sample}.gene.counts.tsv', run=run, sample=SAMPLES_DICT[run]))

## flatten these lists
trimOut = [filename for elem in trimOut for filename in elem]
star_bam_output = [filename for elem in star_bam_output for filename in elem]
verifyBam_output = [filename for elem in verifyBam_output for filename in elem]
picard_RnaSeqMetrics_output = [filename for elem in picard_RnaSeqMetrics_output for filename in elem]
picard_JumpingLibraryMetrics_output = [filename for elem in picard_JumpingLibraryMetrics_output for filename in elem]
star_bam_index_output = [x + '.bai' for x in star_bam_output]
feature_counts_output = [filename for elem in feature_counts_output for filename in elem]
arriba_fusion_gene_output = [filename for elem in arriba_fusion_gene_output for filename in elem]

rule all:
    input:
        trimOut, star_genome_output, star_bam_output, #verifyBam_output, 
        picard_RnaSeqMetrics_output, picard_JumpingLibraryMetrics_output, 
        feature_counts_output, gatk_ase_output, arriba_fusion_gene_output

rule picard_CollectRnaSeqMetrics:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.out.bam'
    output:
        'Data/{run}/star/{sample}/{sample}.Aligned.rna_metrics'
    shell:
        '{picard_cmd} CollectRnaSeqMetrics I={input} O={output} REF_FLAT={refflat_file} STRAND=NONE'

rule picard_CollectJumpingLibraryMetrics:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.jump_metrics'
    shell:
        '{picard_cmd} CollectJumpingLibraryMetrics I={input} O={output}'

rule trim_fastq:
    input:
        fq1="Data/{run}/RAW/{sample}_R1.fastq.gz",
        fq2="Data/{run}/RAW/{sample}_R2.fastq.gz"
    output:
        fq1='Data/{run}/fastq/{sample}_R1_val_1.fq.gz',
        fq2='Data/{run}/fastq/{sample}_R2_val_2.fq.gz'
    priority: 6
    shell:
        '{trim_galore_cmd} --gzip --fastqc '
        '--output_dir Data/{wildcards.run}/fastq '
        '--paired {input.fq1} {input.fq2} '

rule align_with_star_2pass:
    input:
        star_genome_output,
        genome_dir=STAR_GENOME_DIR,
        vcfCommonVariants = dummyVcfStar,
        fq1='Data/{run}/fastq/{sample}_R1_val_1.fq.gz',
        fq2='Data/{run}/fastq/{sample}_R2_val_2.fq.gz'
    output:
        'Data/{run}/star/{sample}/{sample}.Aligned.out.bam'
    params: 
        prefix='Data/{run}/star/{sample}/{sample}.'
    threads: 8
    shell:
        '{star_cmd} --genomeDir {input.genome_dir} '
        '--readFilesIn {input.fq1} {input.fq2} '
        '--outFileNamePrefix {params.prefix} '
        '--outSAMtype BAM Unsorted --waspOutputMode SAMtag ' 
        '--alignSJoverhangMin 8 --varVCFfile {input.vcfCommonVariants} ' 
        '--alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory '
        '--alignIntronMin 20 --alignIntronMax 1000000 '
        '--alignMatesGapMax 1000000 --sjdbScore 2 '
        '--outFilterType BySJout --quantMode TranscriptomeSAM GeneCounts '
        '--outFilterMultimapNmax 20 --outFilterMismatchNmax 999 '
        '--outSAMstrandField intronMotif --chimOutType WithinBAM HardClip '
        '--outFilterIntronMotifs RemoveNoncanonical '
        '--outSAMattributes NH HI NM MD AS XS vA vG vW '
        '--outSAMunmapped Within --twopassMode Basic '
        ' --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 '
        '--alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentMin 10 '
        '--chimScoreDropMax 30 --chimScoreJunctionNonGTAG 0 '
        '--chimScoreSeparation 1 --chimSegmentReadGapMax 3 '
        '--chimMultimapNmax 50 --runThreadN {threads} '
        '--readFilesCommand zcat --chimJunctionOverhangMin 10 ' 

rule sort_alignment:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.out.bam'
    output:
        temp('Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam')
    threads: 4
    priority: 10
    shell:
        'samtools sort -m 12G -@ 4 -O bam -o {output} {input}'

rule write_primary_alignment_bam:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        temp('Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.primary.bam')
    params:
        headerTmp = 'Data/{run}/star/{sample}/header.sam'
    priority: 10
    shell:
        'sh ./getPrimaryAlingments.sh {input} {params.headerTmp} {output}'

rule index_star_bams:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.primary.bam'
    output:
        temp('Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.primary.bam.bai')
    priority: 10
    shell:
        'samtools index {input} '

rule build_star_genome_indexes:
    input:
        fasta=fasta_unzipped,
        annotation=annotation_gtf
    output:
        star_genome_output
    threads: 8
    shell:
        '{star_cmd} --runMode genomeGenerate --genomeDir {STAR_GENOME_DIR} '
        '--genomeFastaFiles {input.fasta} --runThreadN {threads} '
        '--sjdbGTFfile {input.annotation} --sjdbOverhang 74'

rule create_fasta_index:
    input:
        fasta_unzipped
    output:
        fasta_idx
    shell:
        'samtools faidx {input}'

rule create_fasta_dict:
    input:
        fasta_unzipped
    output:
        fasta_dict
    shell:
        '{picard_cmd} CreateSequenceDictionary R={input} O={output}'

rule verify_bam_id:
    '''Verify bam origin, using verify BAM id.'''
    input:
        bam = 'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.primary.bam',
        bai = 'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.primary.bam.bai'
    output:
        'Data/{run}/verifyBamId/{sample}.bestSM'
    priority: 14
    params:
        outBestSm='Data/{run}/verifyBamId/{sample}'
    shell:
        'verifyBamID --vcf /hps/nobackup2/stegle/users/mjbonder/HGSVC/ref/HG38_1kG/ALL.1kG.GRCh38.genotypes.20170504.biallelic.pruned.maf0.1.filtered.vcf.gz '
        '--bam {input.bam} --bai {input.bai} --ignoreRG --out {params.outBestSm} --best --noPhoneHome'

rule count_genes:
    input:
        bam = 'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.bam',
        gtf = annotation_gtf
    output:
        'Data/{run}/initialQuant_featureCounts/{sample}.gene.counts.tsv'
    priority: 4
    run:
        shell("{featurecount_cmd} -B -C -p --primary "
            ' -a {input.gtf} -g gene_id '
            ' -o  {output} {input.bam}')

rule filter_wasp:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam'
    output:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.out.bam'
    threads: 4
    shell:
        'samtools view -h {input} |'
        'python3 /hps/nobackup2/stegle/users/mjbonder/HGSVC/filter_wasp.py |'
        'samtools view -Sb - > {output}'

rule picard_mark_dups:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.out.bam'
    output:
        bam='Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.bam',
        metrics='Data/{run}/star/{sample}/{sample}.output.metrics'
    shell:
        '{picard_cmd} -Xmx50g -Xms50g MarkDuplicates I={input} O={output.bam} CREATE_INDEX=true '
        'VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true M={output.metrics} TMP_DIR=/hps/nobackup2/stegle/tmpFiles/ '

rule picard_read_groups:
    input:
        bam='Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.bam'
    output:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam'
    shell:
        '{picard_cmd}  AddOrReplaceReadGroups I={input.bam} O={output} SO=coordinate '
        'RGID={wildcards.sample} RGLB={wildcards.sample} '
        'RGPL=ILLUMINA RGPU=MACHINE1 RGSM={wildcards.sample}'

rule index_star_ase_bams:
    input:
        'Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam'
    output:
        temp('Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam.bai')
    priority: 10
    shell:
        'samtools index {input} '

rule call_ase_gatk_all:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        dbSnp=mergedVariantCalls,
        bam='Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam',
        bai='Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam.bai'
    output:
        'Data/{run}/ase/{sample}.ase.all.tsv'
    priority: 3
    shell:
        '{gatk_cmd} -T ASEReadCounter -R {input.fasta} -I {input.bam} -o {output} '
        ' -U ALLOW_N_CIGAR_READS -sites {input.dbSnp} '
        '--minMappingQuality 10 --minBaseQuality 2'

rule call_ase_gatk_strelka:
    input:
        fasta=fasta_unzipped,
        fai=fasta_idx,
        fa_dict=fasta_dict,
        dbSnp=strelkaVariants,
        bam='Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam',
        bai='Data/{run}/star/{sample}/{sample}.Aligned.sortedByCoord.WASP.dedup.rgadded.bam.bai'
    output:
        'Data/{run}/ase/{sample}.ase.strelka.tsv'
    priority: 3
    shell:
        '{gatk_cmd} -T ASEReadCounter -R {input.fasta} -I {input.bam} -o {output} '
        ' -U ALLOW_N_CIGAR_READS -sites {input.dbSnp} '
        '--minMappingQuality 10 --minBaseQuality 2'

rule call_fusion_genes_arriba:
    input:
        fasta=fasta_unzipped,
        blacklist=arriba_black_list,
        fa_dict=fasta_dict,
        annotation=annotation_gtf,
        bam='Data/{run}/star/{sample}/{sample}.Aligned.out.bam'
    output:
        'Data/{run}/arriba/{sample}.fusionGenes.tsv'
    priority: 3
    shell:
        '{arriba_cmd} -x {input.bam} -o {output} '
        ' -g {input.annotation} -a {input.fasta} '
        '-b {input.blacklist} '

rule call_fusion_genes_arriba_wgs:
    input:
        fasta=fasta_unzipped,
        blacklist=arriba_black_list,
        fa_dict=fasta_dict,
        annotation=annotation_gtf,
        wgsInfo = arriba_wgs_info,
        bam='Data/{run}/star/{sample}/{sample}.Aligned.out.bam'
    output:
        'Data/{run}/arriba/{sample}.WGS_guided.fusionGenes.tsv'
    priority: 3
    shell:
        '{arriba_cmd} -x {input.bam} -o {output} '
        ' -g {input.annotation} -a {input.fasta} '
        '-b {input.blacklist} -d {input.wgsInfo}'

rule call_fusion_genes_arriba_wgs_noBlackList:
    input:
        fasta=fasta_unzipped,
        fa_dict=fasta_dict,
        annotation=annotation_gtf,
        wgsInfo = arriba_wgs_info,
        bam='Data/{run}/star/{sample}/{sample}.Aligned.out.bam'
    output:
        'Data/{run}/arriba/{sample}.WGS_guided.noBlackList.fusionGenes.tsv'
    priority: 3
    shell:
        '{arriba_cmd} -x {input.bam} -o {output} '
        ' -g {input.annotation} -a {input.fasta} '
        ' -d {input.wgsInfo} -f blacklist'

##arriba -b 

##To add?
#rule split_n_trim_gatk:
#    input:
#        fasta=fasta_unzipped,
#        fai=fasta_idx,
#        fa_dict=fasta_dict,
#        bam='Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.dedup.bam'
#    output:
#        'Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.bam'
#    shell:
#        '{gatk_cmd} -T SplitNCigarReads -R {input.fasta} -I {input.bam} -o {output} '
#        '-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 '
#        '-U ALLOW_N_CIGAR_READS '
#
#rule indel_realignmnet_gatk:
#    input:
#        fasta=fasta_unzipped,
#        bam='Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.bam',
#        targetIntervals=reAlignmentIntervals,
#        known1=knownIndelsMills,
#        known2=knownIndels100G
#    output:
#        'Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bam'
#    shell:
#        '{gatk_cmd} -T IndelRealigner -R {input.fasta} -I {input.bam} '
#        '-targetIntervals {input.targetIntervals} -known {input.known1} -known {input.known2} '
#        '-U ALLOW_N_CIGAR_READS --consensusDeterminationModel KNOWNS_ONLY --LODThresholdForCleaning 0.4  '
#        '-o {output}' 
#
#rule base_recalibrator_gatk:
#    input:
#        fasta=fasta_unzipped,
#        dbSnp=dbSnpVcfBiallelicSnp,
#        bam='Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bam',
#        known1=knownIndelsMills,
#        known2=knownIndels100G
#    output:
#        'Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr'
#    shell:
#        '{gatk_cmd} -T BaseRecalibrator -R {input.fasta} -I {input.bam} '
#        '-knownSites {input.known1} -knownSites {input.known2} -knownSites {input.dbSnp} '
#        '-nct 2 '
#        '-o {output}' 
#
#rule recalibrated_writer_gatk:
#    input:
#        fasta=fasta_unzipped,
#        bam='Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bam',
#        bqsr='Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr'
#    output:
#        'Data/{run}/star/{sample}/{sample}.2pass.Aligned.sortedByCoord.split.realigned.bqsr.bam'
#    shell:
#        '{gatk_cmd} -T PrintReads -R {input.fasta} -I {input.bam} '
#        '-BQSR {input.bqsr} -nct 2 '
#        '-o {output}' 
