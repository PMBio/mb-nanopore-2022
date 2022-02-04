#! /usr/bin/env python

from __future__ import print_function
import argparse
import sys
import cyvcf2
import collections


# Parse command line
parser = argparse.ArgumentParser(description='Phase based on tumor')
parser.add_argument('-v', '--vcf', metavar='input.vcf.gz', required=True, dest='vcf', help='input VCF file (required)')
parser.add_argument('-g', '--germ', metavar='blood', required=True, dest='germ', help='germline sample (required)')
parser.add_argument('-t', '--tumor', metavar='tumor', required=True, dest='tumor', help='tumor sample (required)')
parser.add_argument('-s', '--snpwin', metavar='100', required=True, dest='snpwin', help='#SNPs per window (required)')
parser.add_argument('-o', '--out', metavar='out.vcf', required=True, dest='out', help='output VCF file (required)')
args = parser.parse_args()

# Parameters
snpwin = int(args.snpwin)
afdiff = 0.15
mincov = 20

# Parse VCF
vcf = cyvcf2.VCF(args.vcf, strict_gt=True)
samples = list(vcf.samples)
tumidx = None
if args.tumor:
    if args.tumor not in samples:
        print(args.tumor, "does not exist in VCF!", sep=" ", file=sys.stderr)
        quit()
    tumidx = samples.index(args.tumor)
germidx = None
if args.germ:
    if args.germ not in samples:
        print(args.germ, "does not exist in VCF!", sep=" ", file=sys.stderr)
        quit()
    germidx = samples.index(args.germ)


# Find switch errors
switches = set()
hetpos = set()
if (germidx is not None) and (tumidx is not None):
    rowcount = 0
    afpre = []
    afsuc = []
    diffvec = []
    for record in vcf:
        # Ignore multi-allelics
        if len(record.ALT) > 1:
            continue

        # Ignore non het
        if record.gt_types[germidx] != vcf.HET:
            continue

        # Phased
        gtval = record.genotypes[germidx]
        if not gtval[2]:
            continue

        # Min. coverage check
        ad = record.format('AD')
        if (ad is None) or (ad[tumidx] is None):
            continue
        if sum(ad[tumidx]) <= mincov:
            continue

        # Remember position
        hetpos.add(record.POS)
        
        # Get phased AD
        if gtval[0] == 0: # Ref
            adp = (ad[tumidx][0], ad[tumidx][1])
        else:
            adp = (ad[tumidx][1], ad[tumidx][0])
            
        # Append AF
        if rowcount < snpwin:
            afpre.append((record.POS, adp))
        elif rowcount < 2 * snpwin:
            afsuc.append((record.POS, adp))
        else:
            afpreH1 = 0
            afpreH2 = 0
            afsucH1 = 0
            afsucH2 = 0
            for i in range(snpwin):
                afpreH1 += afpre[i][1][0]
                afpreH2 += afpre[i][1][1]
                afsucH1 += afsuc[i][1][0]
                afsucH2 += afsuc[i][1][1]
            pos = afsuc[rowcount % snpwin][0]
            preval = float(afpreH1) / float(afpreH1 + afpreH2)
            sucval = float(afsucH1) / float(afsucH1 + afsucH2)
            # Switches
            diff = abs(sucval - preval)            
            if (((preval < 0.5) and (sucval > 0.5)) or ((preval > 0.5) and (sucval < 0.5))) and (abs(preval - 0.5) > afdiff / 2.0) and (abs(sucval - 0.5) > afdiff / 2.0) and (diff > afdiff):
                diffvec.append((pos, abs(sucval - preval)))
            else:
                diffvec.append((pos, 0))
            afpre[rowcount % snpwin] = afsuc[rowcount % snpwin]
            afsuc[rowcount % snpwin] = (record.POS, adp)
        rowcount = rowcount + 1

    # Find maxima in diff vector
    localmax = []
    if len(diffvec) > 2:
        for i in range(len(diffvec)-1):
            if i > 0:
                if (diffvec[i][1] > diffvec[i-1][1]) and (diffvec[i][1] > diffvec[i+1][1]):
                    localmax.append(diffvec[i])

    # Find best maxima in SNP * 1000bp distance
    offset = snpwin * 1000
    for i in range(len(localmax)):
        bettermax = False
        for j in range(len(localmax)):
            if i != j:
                if abs(localmax[i][0] - localmax[j][0]) <= offset:
                    if (localmax[j][1] > localmax[i][1]) or ((localmax[j][1] == localmax[i][1]) and (j > i)):
                        bettermax = True
                        break
        if not bettermax:
            print(localmax[i])
            switches.add(localmax[i][0])
vcf.close()

# Switch phase at switch points
print("#Switches: ", len(switches))
print(sorted(switches))
if True:
    vcfIn = cyvcf2.VCF(args.vcf, strict_gt=True)
    w = cyvcf2.Writer(args.out, vcfIn)
    switch = False
    for record in vcfIn:
        if record.POS not in hetpos:
            continue

        if record.POS in switches:
            switch = not switch

        # Switch if necessary
        if switch:
            gtlist = record.genotypes
            gtval = record.genotypes[germidx][0]
            gtlist[germidx][0] = gtlist[germidx][1]
            gtlist[germidx][1] = gtval
            record.genotypes = gtlist

        # Keep site
        w.write_record(record)
            
    # Close file
    w.close()
    vcfIn.close()
