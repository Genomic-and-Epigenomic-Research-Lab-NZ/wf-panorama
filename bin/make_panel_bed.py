#!/usr/bin/env python3
"""
CLI wrapper for make_panel_bed functionality adapted from workflow/scripts/make_panel_bed.py
"""
import argparse
import os
import pandas as pd
import numpy as np

from warnings import WarningMessage


def add_functional_flanking_regions(panel_csv):
    promoter_length = 2000
    end_length = 1000

    for index, entry in panel_csv.iterrows():
        if entry['Is this record the whole gene?'] == 'Yes':
            if isinstance(entry['start pos'], str) and ',' in entry['start pos']:
                start_coord = int(entry['start pos'].replace(',', ''))
            else:
                start_coord = int(entry['start pos'])
            if isinstance(entry['end pos'], str) and ',' in entry['end pos']:
                end_coord = int(entry['end pos'].replace(',', ''))
            else:
                end_coord = int(entry['end pos'])

            new_start_coord = None
            new_end_coord = None
            strand = entry['strand']

            if strand == '+':
                gene_start_coord = start_coord
                gene_end_coord = end_coord
                new_start_coord = gene_start_coord - promoter_length
                new_end_coord = gene_end_coord + end_length
            elif strand == '-':
                gene_start_coord = end_coord
                gene_end_coord = start_coord
                new_start_coord = gene_end_coord - end_length
                new_end_coord = gene_start_coord + promoter_length

            panel_csv.at[index, 'start pos'] = new_start_coord
            panel_csv.at[index, 'end pos'] = new_end_coord

        if entry['strand'] not in ["+", "-"]:
            panel_csv.at[index, 'strand'] = "*"

    return panel_csv


def restructure_to_bed(df):
    df['score'] = np.nan
    df = df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]
    df = df.copy()
    df['chromStart'] = df['chromStart'].round().astype(int)
    df['chromEnd'] = df['chromEnd'].round().astype(int)
    return df


def get_immune_infiltrate_ref_locs(immune_reference_dataset, epic_locs_hg38):
    imm_ref = pd.read_csv(immune_reference_dataset)
    if 'NAME' in imm_ref.columns:
        imm_ref.rename(columns={'NAME': 'CpGs'}, inplace = True)
    probe_list = imm_ref['CpGs']

    EPIC_probes_loc_hg38 = pd.read_csv(epic_locs_hg38, index_col=0)
    probes_genomic_locs = EPIC_probes_loc_hg38[EPIC_probes_loc_hg38['probe'].isin(probe_list)]

    if len(probes_genomic_locs) < len(probe_list):
        WarningMessage(
            f'not all genomic locations were found for reference probes ({len(probe_list) - len(probes_genomic_locs)} are missing.)',
            category = UserWarning,
            filename = 'make_panel_bed.py',
            lineno = 71)

    immune_probes_bed_df = probes_genomic_locs.rename(columns = {
        'probe': 'name',
        'seqnames': '#chrom',
        'start': 'chromStart',
        'end': 'chromEnd'
    })
    immune_probes_bed_df = restructure_to_bed(immune_probes_bed_df)

    return immune_probes_bed_df


def add_immune_infiltrate_locations(input_bed, immune_reference_dataset, epic_locs_hg38):
    immune_probes_bed_df = get_immune_infiltrate_ref_locs(immune_reference_dataset, epic_locs_hg38)
    merged_df = pd.concat([input_bed, immune_probes_bed_df])
    return merged_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create panel bed and all targets from panel_metadata.csv')
    parser.add_argument('--panel-csv', required=True)
    parser.add_argument('--immune-reference', required=True)
    parser.add_argument('--epic-locs', required=True)
    parser.add_argument('--panel-bed', required=True)
    parser.add_argument('--all-targets', required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.panel_bed), exist_ok=True)
    os.makedirs(os.path.dirname(args.all_targets), exist_ok=True)

    panel_csv = pd.read_csv(args.panel_csv, dtype={'ID': str}, thousands = ',')
    panel_csv = panel_csv[~panel_csv.ID.str.startswith("4")]
    panel_csv = panel_csv[~panel_csv.ID.str.startswith("5")]

    panel_csv = add_functional_flanking_regions(panel_csv)

    panel_bed = panel_csv[['ID', 'chrom', 'start pos', 'end pos', 'strand']]
    panel_bed = panel_bed.rename(columns = {
        'ID': 'name', 
        'chrom': '#chrom', 
        'start pos': 'chromStart', 
        'end pos':'chromEnd'})
    panel_bed = restructure_to_bed(panel_bed)

    panel_bed.to_csv(args.panel_bed, sep = '\t', index = False, header = False)

    all_targets = add_immune_infiltrate_locations(panel_bed, args.immune_reference, args.epic_locs)
    all_targets.to_csv(args.all_targets, sep = '\t', index = False, header = False)
