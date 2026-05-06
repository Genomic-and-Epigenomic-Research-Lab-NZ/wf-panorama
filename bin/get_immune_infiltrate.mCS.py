#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np
from shared_functions import variant_prep, BIOMARKER_ID, BIOMARKER_NAME, BIOMARKER_TYPE, BIOMARKER_TYPE_FULL, SCORING_TYPE, RESULT_OPTIONS, preclin_stage_panel_result_header

lymphocytes = ["CD19", "CD4_Eff",
                "CD56", "CD8", "Treg", ]
monocytes = ["CD14"]
neutrophils = ["Neu"]
eosinophils = ["Eos"]


def get_LMR(deconv_df):
    """Lymphocyte/Monocyte ratio"""
    lymphocyte_val = sum(deconv_df[lymphocytes].iloc[0])
    monocyte_val = sum(deconv_df[monocytes].iloc[0])
    # Override for testing
    # monocyte_val = 0.015
    if monocyte_val < 0.0001:
        print("Monocyte value too low to get LMR.", file=sys.stderr)
        LMR_ratio = np.nan
    elif lymphocyte_val < 0.0001:
        print("Lymphocyte value too low to get LMR.", file=sys.stderr)
        LMR_ratio = np.nan
    else:
        LMR_ratio = lymphocyte_val / monocyte_val
        print(f"Lymphocyte to Monocyte ratio (LMR) is {LMR_ratio:.2f}", file=sys.stderr)
    return LMR_ratio


def get_NLR(deconv_df):
    """Neutrophil/Lymphocyte ratio"""
    lymphocyte_val = sum(deconv_df[lymphocytes].iloc[0])
    neutrophil_val = sum(deconv_df[neutrophils].iloc[0])
    if neutrophil_val < 0.0001:
        print("Neutrophil value too low to get NLR.", file=sys.stderr)
        NLR_ratio = np.nan
    elif lymphocyte_val < 0.0001:
        print("Lymphocyte value too low to get NLR.", file=sys.stderr)
        NLR_ratio = np.nan
    else:
        NLR_ratio = neutrophil_val / lymphocyte_val
        print(f"Neutrophil to Lymphocyte ratio (NLR) is {NLR_ratio:.2f}", file=sys.stderr)
    return NLR_ratio


def get_dendritic_cells():
    """Future"""
    return


def get_Th1_cells():
    """Future"""
    return


def get_M2_macrophages():
    """Future"""
    return


def get_platelet_count():
    """Future"""
    return


parser = argparse.ArgumentParser()
parser.add_argument('--deconv', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--panel', required=True)
args = parser.parse_args()

deconv_df = pd.read_csv(args.deconv, sep=',')  # TODO check once CS is up

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

# Get ratios of interest
LMR_ratio = get_LMR(deconv_df)
NLR_ratio = get_NLR(deconv_df)

# Placeholders
# LMR_ratio = np.nan
# NLR_ratio = np.nan

panel_data_ratio = variant_prep(args.panel, 'immune_ratio')
panel_data_infiltrate = variant_prep(args.panel, 'immune_inf')

bm_classif_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)

for i, row in panel_data_ratio.iterrows():
    if row[BIOMARKER_NAME] == 'LMR':
        result = LMR_ratio
    elif row[BIOMARKER_NAME] == 'NLR':
        result = NLR_ratio
    else:
        result = np.nan
    bm_classif_panel_df.loc[i] = [row[BIOMARKER_ID], row[BIOMARKER_NAME], row[SCORING_TYPE], row[BIOMARKER_TYPE_FULL], row[RESULT_OPTIONS], result]

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
    
for i, row in panel_data_infiltrate.iterrows():
    name = row[BIOMARKER_NAME]
    result = mapping.get(name, np.nan)
    bm_classif_panel_df.loc[i] = [row[BIOMARKER_ID], row[BIOMARKER_NAME], row[SCORING_TYPE], row[BIOMARKER_TYPE_FULL], row[RESULT_OPTIONS], result]

bm_classif_panel_df.to_csv(args.out, index=False)
