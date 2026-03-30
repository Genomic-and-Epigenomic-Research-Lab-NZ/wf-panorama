#!/usr/bin/env python3
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True)
parser.add_argument('--output', required=True)
args = parser.parse_args()

header = ["chrom", "start", "end", "mod base code", "score", "strand", "start pos", "end pos", "color", "N valid cov", "fraction modified", "N mod", "N canonical", "N other mod", "N delete", "N fail", "N diff", "N nocall"]

with open(args.input, 'r') as fb:
    bed_df = pd.read_csv(fb, sep='\t', names=header)

dss_df = bed_df.copy()[['chrom', 'start', 'N valid cov', 'N mod']]
dss_df.rename({'chrom':'chr', 'start': 'pos', 'N valid cov': 'N', 'N mod': 'X'}, axis=1, inplace=True)

dss_df.to_csv(args.output, sep='\t', index=False)
