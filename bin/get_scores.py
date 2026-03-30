#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
from shared_functions import BIOMARKER_TYPE

parser = argparse.ArgumentParser()
parser.add_argument('--panel', required=True)
parser.add_argument('--snv', required=True)
parser.add_argument('--mod', required=True)
parser.add_argument('--immune', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()

panel_df = pd.read_csv(args.panel, dtype={"ID": str})
with open(args.snv,'r') as f: snv_df = pd.read_csv(f, dtype={"ID": str})
with open(args.mod,'r') as f: mod_df = pd.read_csv(f, dtype={"ID": str})
with open(args.immune,'r') as f: immune_df = pd.read_csv(f, dtype={"ID": str})

panel_df = panel_df[panel_df['Is there either data?'] == 'Yes']

results_dict = {}
min_scores_list = []
max_scores_list = []

for _, panel_entry in panel_df.iterrows():
    panel_id = panel_entry['ID']
    if panel_id in snv_df['ID'].values:
        result_entry = snv_df[snv_df['ID']==panel_id].iloc[0]
    elif panel_id in mod_df['ID'].values:
        result_entry = mod_df[mod_df['ID']==panel_id].iloc[0]
    elif panel_id in immune_df['ID'].values:
        result_entry = immune_df[immune_df['ID']==panel_id].iloc[0]
    else:
        continue
    # Simplified scoring: map result to 0-100 depending on type
    if panel_entry[BIOMARKER_TYPE] == 'immune_inf':
        val = float(result_entry.get('Result',0))/100.0
    else:
        val = float(result_entry.get('Result',0))
    score = val * 100
    results_dict[panel_id] = [result_entry.get('Result', ''), np.float64(score)]

# Aggregate
total_score = np.nansum([v[1] for v in results_dict.values()])
results_dict['Total (raw)'] = ['', total_score]
# normalise placeholder
results_dict['Total (normalised to 0.0 - 100.0)'] = ['', (total_score/100.0)*100 if total_score else 0]

score_results_df = pd.DataFrame.from_dict(results_dict, orient='index', columns=['Result','Score'])
score_results_df.to_csv(args.out)
