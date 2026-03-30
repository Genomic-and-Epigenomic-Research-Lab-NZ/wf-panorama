#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
from shared_functions import variant_prep, BIOMARKER_NAME, SCORING_TYPE, RESULT_OPTIONS

parser = argparse.ArgumentParser()
parser.add_argument('--deconv', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--panel', required=True)
args = parser.parse_args()

deconv_df = pd.read_csv(args.deconv, sep='\t')

Monocytes = deconv_df['CD14'].iloc[0]
Bcells = deconv_df['CD19'].iloc[0]
CD4_Tcells = deconv_df['CD4_Eff'].iloc[0]
NK_cells = deconv_df['CD56'].iloc[0]
CD8_Tcells = deconv_df['CD8'].iloc[0]
Tregs = deconv_df['Treg'].iloc[0]
Endothelial = deconv_df['Endothelial'].iloc[0]
Eosinophils = deconv_df['Eos'].iloc[0]
Fibroblasts = deconv_df['Fibroblast'].iloc[0]
Neutrophils = deconv_df['Neu'].iloc[0]
Cancer = deconv_df['Cancer'].iloc[0]

LMR_ratio = np.nan
NLR_ratio = np.nan

panel_data_ratio = variant_prep(args.panel, 'immune_ratio')
panel_data_infiltrate = variant_prep(args.panel, 'immune_inf')

bm_classif_panel_df = pd.DataFrame(columns=['ID', 'Biomarker name', 'Scoring Type', 'Biomarker Type', 'Result Options', 'Result'])

for i, row in panel_data_ratio.iterrows():
    if row[BIOMARKER_NAME] == 'LMR':
        result = LMR_ratio
    elif row[BIOMARKER_NAME] == 'NLR':
        result = NLR_ratio
    else:
        result = np.nan
    bm_classif_panel_df.loc[i] = [row['ID'], row[BIOMARKER_NAME], row[SCORING_TYPE], row['Biomarker Type'], row[RESULT_OPTIONS], result]

for i, row in panel_data_infiltrate.iterrows():
    name = row[BIOMARKER_NAME]
    mapping = {
        'Monocyte_inf': Monocytes,
        'Bcell_inf': Bcells,
        'CD4_inf': CD4_Tcells,
        'NK_inf': NK_cells,
        'CD8_inf': CD8_Tcells,
        'Treg_inf': Tregs,
        'Neutrophil_inf': Neutrophils,
        'Endothelial_inf': Endothelial,
        'Eosinophil_inf': Eosinophils,
        'Fibroblast_inf': Fibroblasts,
        'Cancer_inf': Cancer
    }
    result = mapping.get(name, np.nan)
    bm_classif_panel_df.loc[i] = [row['ID'], row[BIOMARKER_NAME], row[SCORING_TYPE], row['Biomarker Type'], row[RESULT_OPTIONS], result]

bm_classif_panel_df.to_csv(args.out, index=False)
