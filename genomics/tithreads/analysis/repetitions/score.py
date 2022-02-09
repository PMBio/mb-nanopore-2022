#! /usr/bin/env python

from __future__ import print_function
from readfq import readfq
import csv
import argparse
import collections
import os
import json
import gzip

# Parse command line
parser = argparse.ArgumentParser(description='Alignment parser')
parser.add_argument('-a', '--alignment', required=True, dest='align', help='templated insertions (required)')
parser.add_argument('-m', '--match', required=True, dest='match', help='match (required)')
parser.add_argument('-n', '--mismatch', required=True, dest='mismatch', help='mismatch (required)')
parser.add_argument('-g', '--gapopen', required=True, dest='gapopen', help='gapopen (required)')
parser.add_argument('-e', '--gapextend', required=True, dest='gapextend', help='gapextend (required)')
args = parser.parse_args()

match = int(args.match)
mismatch = int(args.mismatch)
gep = int(args.gapextend)
gop = int(args.gapopen)

with (gzip.open(args.align) if args.align.endswith('.gz') else open(args.align)) as f:
    h1 = None
    h2 = None
    for seqName, seqNuc, seqQuals in readfq(f):
        if h1 is None:
            h1 = seqName
            h1 = h1.replace('>','')
            f1 = seqNuc
        else:
            h2 = seqName
            h2 = h2.replace('>','')
            f2 = seqNuc
    # Get alignment start and end
    leadGap = True
    seqlen = 0
    seq1start = 0
    for idx, (c1, c2) in enumerate(zip(f1, f2)):
        if (c2 == '-'):
            pass
        else:
            seqlen += 1
            alnEnd = idx
            if leadGap:
                alnStart = idx
                leadGap = False
        if c1 == '-':
            pass
        else:
            if leadGap:
                seq1start += 1

                
    # Score alignment
    score = 0
    ingap = False;
    seq1end = seq1start
    for idx, (c1, c2) in enumerate(zip(f1, f2)):
        if (c2 == '-'):
            if (idx >= alnStart) and (idx <= alnEnd):
                if ingap:
                    score += gep
                else:
                    ingap = True
                    score += gop + gep
        else:
            if c1 == '-':
                if (idx >= alnStart) and (idx <=alnEnd):
                    if ingap:
                        score += gep
                    else:
                        ingap = True
                        score += gop + gep
            else:
                ingap = False
                if (idx >= alnStart) and (idx <= alnEnd):
                    seq1end += 1
                    if (c1 == c2):
                        score += match
                    else:
                        score += mismatch

# Length-normalized score
if score < 0:
    score = -1
else:
    score = int(float(score) / float(seqlen * match) * 100)
    if score < 40:  # <80% identity threshold
        score = -1
print(h1, seq1start, seq1end, score, sep="\t")
