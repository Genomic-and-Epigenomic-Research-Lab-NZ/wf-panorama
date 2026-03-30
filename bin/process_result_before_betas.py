#!/usr/bin/env python3
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser()
parser.add_argument('--epic', required=True)
parser.add_argument('--dss', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()

EPIClocs_df = pd.read_csv(args.epic, sep=',')
dss_df = pd.read_csv(args.dss, sep='\t')

def filter_dss_to_x_depth(dss_df, depth = 30):
    dss_df = dss_df[dss_df['N'] >= depth]
    dss_df['endpos'] = dss_df['pos'] + 1
    return dss_df


def add_illumina_probes(dss_df, EPIClocs_df):
    merged_df = pd.merge(EPIClocs_df, dss_df, left_on=['seqnames','start'], right_on=['chr','pos'], how='right')
    merged_df_clean = merged_df[merged_df['N'].notna()].copy()
    merged_df_clean['pos'] = merged_df_clean['pos'].astype(int)
    merged_df_clean['N'] = merged_df_clean['N'].astype(int)
    merged_df_clean['X'] = merged_df_clean['X'].astype(int)
    merged_df_clean.drop('start', axis=1, inplace=True)
    df = merged_df_clean[['probe', 'chr','pos','strand','N','X']]
    return df

os.makedirs(os.path.dirname(args.out), exist_ok=True)

# Use depth=1 for development as original had TODO
dss_df_30x = filter_dss_to_x_depth(dss_df.copy(), depth=1)

dss_df_30x_probes = add_illumina_probes(dss_df_30x, EPIClocs_df)

dss_df_30x_probes.to_csv(args.out, index=False)
