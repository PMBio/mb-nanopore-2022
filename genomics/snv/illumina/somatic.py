#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import sys
import collections
import json
import cyvcf2

# Parse command line
parser = argparse.ArgumentParser(description='AF comparison tumor vs. relapse')
parser.add_argument('-v', '--variants', metavar='sample.vcf.gz', required=True, dest='variants', help='VCF file (required)')
args = parser.parse_args()

min_cov = 20
min_af = 0.1
vaf = dict()

# Desired VEP columns
descols = ["Consequence", "IMPACT", "SYMBOL", "Gene", "STRAND", "Feature", "EXON", "HGVSc", "HGVSp", "cDNA_position", "CDS_position", "Protein_position", "Amino_acids", "Codons", "Existing_variation", "SIFT", "PolyPhen", "gnomADg_AF", "ExAC_AF", "ExAC_nonTCGA_AF", "gnomAD_genomes_non_cancer_AF", "MAX_AF", "CLIN_SIG", "CADD_phred_hg19", "phastCons100way_vertebrate", "phastCons17way_primate"]
selected = ['missense_variant', 'splice_acceptor_variant', 'splice_donor_variant', 'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant', 'inframe_insertion', 'inframe_deletion', 'synonymous_variant']

# VCF/BCF parsing
varstore = collections.defaultdict(collections.defaultdict)
vcf = cyvcf2.VCF(args.variants)
for header in vcf.header_iter():
    if header['HeaderType'] == "INFO":
        if header['ID'] == "CSQ":
            vepcols = header['Description'].replace('Consequence annotations from Ensembl VEP. Format: "', '')[:-1].split('|')
            #print('\n'.join(vepcols))

addr = dict()
for idx, val in enumerate(vepcols):
    addr[val] = idx

print("chromosome", "position", "REF", "ALT", "tumor_vaf", "relapse_vaf", '\t'.join(descols), sep='\t')
#print("chromosome", "position", "REF", "ALT", "tumor_vaf", "relapse_vaf", sep='\t')
for record in vcf:
    if len(record.ALT) > 1:
        continue
    ad = record.format('AD')
    dp = record.format('DP')
    for idx, spl in enumerate(vcf.samples):
        vaf[spl] = 0
        if dp[idx] >= min_cov:
            vaf[spl] = float(ad[idx][1])/float(ad[idx][0] + ad[idx][1])
    #if (vaf['tumor'] > min_af) or (vaf['relapse'] > min_af):
        #print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), vaf['tumor'], vaf['relapse'], sep='\t')
    #continue
    csq = record.INFO.get('CSQ')
    if csq is None:
        print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), "unknown_variant", sep='\t')
        continue
    transcripts = csq.split(',')
    for tr in transcripts:
        fields = tr.split('|')
        consout = False
        for constype in fields[addr['Consequence']].split('&'):
            if constype in selected:
                consout = True
                break
        af_vcf = record.INFO.get('AF')
        if (consout) and (fields[addr['CANONICAL']] == "YES") and (fields[addr['BIOTYPE']] == "protein_coding"):
            if len(fields) != len(vepcols):
                print("Error: VEP annotation is corrupted!", file=sys.stderr)
                sys.exit(1)
            fields[addr['EXON']] = fields[addr['EXON']].replace('/', ';')
            if (vaf['tumor'] > min_af) or (vaf['relapse'] > min_af):
                print(record.CHROM, record.POS, record.REF, ','.join(record.ALT), vaf['tumor'], vaf['relapse'], '\t'.join([fields[addr[cname]] for cname in descols]), sep='\t')

