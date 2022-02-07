#! /usr/bin/env python

from __future__ import print_function
import csv
import argparse
import collections
import os
import json
import gzip

# Parse command line
parser = argparse.ArgumentParser(description='Parse MUMmer matches')
parser.add_argument('-m', '--matches', required=True, dest='matches', help='MUMmer matches (required)')
args = parser.parse_args()

print("ref", "xstart", "xend", "query", "ystart", "yend", "direction", sep="\t")
if (os.path.exists(args.matches)) and (os.path.isfile(args.matches)):
    with open(args.matches) as f:
        header = None
        reverse = False
        for line in f:
            content = line.strip()
            if content.startswith(">"):
                header = content.split()
                header = header[1]
                if content.endswith("Reverse"):
                    reverse = True
                else:
                    reverse = False
            else:
                fields = content.split()
                rpos, qpos, l = [int(x) for x in fields[1:]]
                if reverse:
                    print(fields[0], rpos, rpos+l-1, header, qpos, qpos-l+1, 'reverse', sep="\t")
                else:
                    print(fields[0], rpos, rpos+l-1, header, qpos, qpos+l-1, 'forward', sep="\t")
