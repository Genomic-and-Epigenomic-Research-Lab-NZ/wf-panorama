#!/usr/bin/env python3
import argparse
import pandas as pd
from shared_functions import variant_prep, preclin_stage_panel_result_header

parser = argparse.ArgumentParser()
parser.add_argument('--panel', required=True)
parser.add_argument('--snv', required=True)
parser.add_argument('--mod', required=True)
parser.add_argument('--immune', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()

snv_df = pd.read_csv(args.snv, dtype={"ID": str})
mod_df = pd.read_csv(args.mod, dtype={"ID": str})
immune_df = pd.read_csv(args.immune, dtype={"ID": str})

preclin_stage_panel_result_header = ["ID", "Marker name", "Scoring Type", "Biomarker Type", "Result Options", "Result"]
full_panel_results_df = pd.concat([snv_df, mod_df, immune_df])
full_panel_results_df = full_panel_results_df.sort_values(['ID'])
full_panel_results_df.reset_index(drop=True, inplace=True)

full_panel_results_df.to_csv(args.out, index=False)
