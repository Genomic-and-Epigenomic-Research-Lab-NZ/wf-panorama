#!/usr/bin/env python3
"""
CLI wrapper for make_panel_bed functionality adapted from workflow/scripts/make_panel_bed.py
"""
import argparse
import pandas as pd
import numpy as np
from collections import defaultdict

from warnings import WarningMessage


def parse_panel(panel_csv_fp):
    """
    Parse the panel CSV file and return a DataFrame with the relevant columns.
    """
    
    panel_dtypes = defaultdict(lambda: "string", {
        "start pos": np.float64,
        "end pos": np.float64,
        "length": np.float64,
    })

    panel_csv = pd.read_csv(panel_csv_fp, dtype=panel_dtypes, thousands = ',')
    
    # Check for decimals in coordinate columns
    # coord_cols = ['start pos', 'end pos']
    # for col in coord_cols:
    #     if col in panel_csv.columns:
    #         has_decimals = panel_csv[col].astype(str).str.contains(r'\.', regex=True).any()
    #         if has_decimals:
    #             raise ValueError(
    #                 f"Error: Column '{col}' contains decimal values.\n" +
    #                 f"Coordinate columns must be integers (whole numbers only).\n" +
    #                 f"Please check your input file ({args.panel_csv}) and ensure all coordinate values are whole numbers."
    #             )
    
    panel_csv = panel_csv[~panel_csv.ID.str.startswith("3")]  # combined markers
    panel_csv = panel_csv[~panel_csv.ID.str.startswith("4")]  # immune infiltrate markers
    panel_csv = panel_csv[~panel_csv.ID.str.startswith("5")]  # non-molecular markers
    
    panel_csv = panel_csv.astype({"start pos": "int64", "end pos": "int64", "length": "int64"})
    
    return panel_csv


def add_msk_impact(panel_csv, msk_impact):
    """msk_impact is a bed file"""
    
    msk_panel = pd.read_csv(msk_impact, sep="\t", header=None, names=["chrom", "start pos", "end pos", "ID", "score", "strand"])
    msk_panel["Is this record the whole gene?"] = "Yes"
    msk_panel["length"] = msk_panel["end pos"] - msk_panel["start pos"]
    
    # add the msk panel to the existing panel_csv
    # note that need to match the columns to the order in panel_csv, they are in different order in msk_panel
    msk_panel = msk_panel[["ID", "chrom", "start pos", "end pos", "strand", "Is this record the whole gene?", "length"]]
    panel_csv = pd.concat([panel_csv, msk_panel], ignore_index=True)
    
    return panel_csv


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


def restructure_to_bed(df, panel = False):
    
    if panel == True:
        df = df[['ID', 'chrom', 'start pos', 'end pos', 'strand']]
        df = df.rename(columns = {
            'ID': 'name', 
            'chrom': '#chrom', 
            'start pos': 'chromStart', 
            'end pos':'chromEnd'})
        
    df['score'] = np.nan
    df = df[['#chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']]
    # df = df.copy()
    # df['chromStart'] = df['chromStart'].round().astype(int)
    # df['chromEnd'] = df['chromEnd'].round().astype(int)
    
    return df


def add_immune_infiltrate_locations(input_bed, immune_reference_dataset, epic_locs_hg38, allEPIC):
    
    EPIC_probes_loc_hg38 = pd.read_csv(epic_locs_hg38, index_col=0)
    # header = "probe","seqnames","start","end","width","strand"
    # row = "cg18478105","chr20",63216298,63216298,1,"-"
    
    if allEPIC:
        probes_bed_df = EPIC_probes_loc_hg38.rename(
            columns={
                'probe': 'name',
                'seqnames': '#chrom',
                'start': 'chromStart',
                'end': 'chromEnd',
                'width': 'length'
                }
            )
    else:
        # Override the immune_reference_dataset found by cancer type
        # immune_reference_dataset = "/path/to/manual/override/immune_reference.csv"
        
        imm_ref = pd.read_csv(immune_reference_dataset)
        if 'NAME' in imm_ref.columns:
            imm_ref.rename(columns={'NAME': 'CpGs'}, inplace = True)
        probe_list = imm_ref['CpGs']

        probes_genomic_locs = EPIC_probes_loc_hg38[EPIC_probes_loc_hg38['probe'].isin(probe_list)]

        if len(probes_genomic_locs) < len(probe_list):
            WarningMessage(
                f'not all genomic locations were found for reference probes ({len(probe_list) - len(probes_genomic_locs)} are missing.)',
                category = UserWarning,
                filename = 'make_panel_bed.py',
                lineno = 71)

        probes_bed_df = probes_genomic_locs.rename(
            columns={
                'probe': 'name',
                'seqnames': '#chrom',
                'start': 'chromStart',
                'end': 'chromEnd',
                'width': 'length'
                }
            )
        
        probes_bed_df = restructure_to_bed(probes_bed_df)        
    
    merged_df = pd.concat([input_bed, probes_bed_df])
    
    # delete length column
    if 'length' in merged_df.columns:
        merged_df = merged_df.drop(columns=['length'])
        
    return merged_df


def main(args):
    
    panel_csv = parse_panel(args.panel_csv)
    
    if args.msk_impact:
        panel_csv = add_msk_impact(panel_csv, args.msk_impact)

    panel_csv = add_functional_flanking_regions(panel_csv)
    
    panel_bed = restructure_to_bed(panel_csv, panel=True)

    panel_bed.to_csv(args.panel_bed, sep = '\t', index = False, header = False)

    all_targets = add_immune_infiltrate_locations(panel_bed, args.immune_reference, args.epic_locs, args.allEPIC)
    
    all_targets.to_csv(args.all_targets, sep = '\t', index = False, header = False)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create panel bed and all targets from panel_metadata.csv')
    parser.add_argument('--panel-csv', required=True)
    parser.add_argument('--immune-reference', required=True)
    parser.add_argument('--epic-locs', required=True)
    parser.add_argument('--panel-bed', required=True)
    parser.add_argument('--all-targets', required=True)
    parser.add_argument('--all-epic', '--allEPIC', dest='allEPIC', action='store_true', help='If set, will add all EPIC probes to the all targets bed file instead of just the immune infiltrate reference probes. WARNING: This will make the targets file cover ~20% of the genome when using a buffersize of 10kb')
    parser.add_argument('--msk_impact', help='If given, will add the msk_impact panel to the panel bed file. WARNING: This will make the panel bed file cover ~70% of the genome when using a buffersize of 10kb')
    args = parser.parse_args()

    main(args)
