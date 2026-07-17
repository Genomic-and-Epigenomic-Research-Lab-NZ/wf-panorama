#!/usr/bin/env python3
import argparse
import pandas as pd
import os
import numpy as np
from shared_functions import preclin_stage_panel_result_header, variant_prep, CHROMOSOMES, SCORING_TYPE, RESULT_OPTIONS, BIOMARKER_NAME, get_BM_TYPE_FULL

# Note: this is a conversion of the original modification_calling.py adapted to CLI usage.

parser = argparse.ArgumentParser()
parser.add_argument('--panel', required=True)
parser.add_argument('--mod_data', required=True)
parser.add_argument('--out-meth', required=True)
parser.add_argument('--out-mod', required=True)
parser.add_argument('--out-raw', required=True)
args = parser.parse_args()

# Basic helper functions reproduced/adapted from original script

def export_for_methatlas(df, fp):
    df = df[df['probe'].notna()]
    df = df.drop(columns=[c for c in ['chr','pos','strand','N','X'] if c in df.columns], errors='ignore')
    df.to_csv(fp, index=False)


def add_prom_start_end(panel_metadata_df_meth):
    def calc_prom_start(row):
        if row.get('Is this record the whole gene?') == 'Yes':
            if row.get('strand') == '+':
                return int(str(row['start pos']).replace(',','')) - 2000
            else:
                return int(str(row['end pos']).replace(',','')) - 1000
        else:
            return int(str(row['start pos']).replace(',','')) - 5
    def calc_prom_end(row):
        if row.get('Is this record the whole gene?') == 'Yes':
            if row.get('strand') == '+':
                return int(str(row['end pos']).replace(',','')) + 1000
            else:
                return int(str(row['start pos']).replace(',','')) + 2000
        else:
            return int(str(row['end pos']).replace(',','')) + 5
    panel_metadata_df_meth['prom_start'] = panel_metadata_df_meth.apply(calc_prom_start, axis=1)
    panel_metadata_df_meth['prom_end'] = panel_metadata_df_meth.apply(calc_prom_end, axis=1)
    return panel_metadata_df_meth


def add_downstream_start_end(panel_metadata_df_meth):
    def calc_down_start(row):
        if row.get('Is this record the whole gene?') == 'Yes':
            if row.get('strand') == '+':
                return int(str(row['start pos']).replace(',','')) - 1000
            else:
                return int(str(row['end pos']).replace(',','')) - 5
        else:
            return int(str(row['start pos']).replace(',','')) - 5
    def calc_down_end(row):
        if row.get('Is this record the whole gene?') == 'Yes':
            if row.get('strand') == '+':
                return int(str(row['end pos']).replace(',','')) + 2000
            else:
                return int(str(row['start pos']).replace(',','')) + 1000
        else:
            return int(str(row['end pos']).replace(',','')) + 5
    panel_metadata_df_meth['down_start'] = panel_metadata_df_meth.apply(calc_down_start, axis=1)
    panel_metadata_df_meth['down_end'] = panel_metadata_df_meth.apply(calc_down_end, axis=1)
    return panel_metadata_df_meth


def sort_and_prep_dfs(variants_df, dss_betas_df):
    # set chromosome column types
    variants_df.chrom = variants_df.chrom.astype("category")
    variants_df.chrom = variants_df.chrom.cat.set_categories(CHROMOSOMES)
    dss_betas_df.chr = dss_betas_df.chr.astype("category")
    dss_betas_df.chr = dss_betas_df.chr.cat.set_categories(CHROMOSOMES)

    # convert locations to ints
    variants_df["prom_start"] = variants_df["prom_start"].astype("int64")
    variants_df["prom_end"] = variants_df["prom_end"].astype("int64")
    variants_df["down_start"] = variants_df["down_start"].astype("int64")
    variants_df["down_end"] = variants_df["down_end"].astype("int64")
    variants_df["start pos"] = (
        variants_df["start pos"].str.replace(",", "").astype("int64")
    )
    variants_df["end pos"] = (
        variants_df["end pos"].str.replace(",", "").astype("int64")
    )
    dss_betas_df["pos"] = dss_betas_df["pos"].astype("int64")

    # Sort by columns: 'chrom' (ascending), 'start pos' (ascending), 'end pos' (ascending)
    variants_df = variants_df.sort_values(["chrom", "prom_start", "prom_end"])
    dss_betas_df = dss_betas_df.sort_values(["chr", "pos"])
    return variants_df, dss_betas_df


def merge_data_to_list_dicts(results_df, panel_data):
    results_as_list_dicts = list()
    for group_chr in results_df.groupby(by="chr", observed=False):
        chrom = group_chr[0]
        # log.write(f"Processing chromosome: {chrom}\n")
        # if chrom not in set(panel_data["chrom"]):
        #     log.write("    Chromosome not in panel regions of interest; skipping.\n")
        #     continue

        # Get chromosome group from variants
        df = panel_data.copy()
        chr_entries = df.loc[df["chrom"] == group_chr[0]].copy()
        chr_entries["min"] = chr_entries[
            ["prom_start", "prom_end", "down_start", "down_end"]
        ].min(axis=1)
        chr_entries["max"] = chr_entries[
            ["prom_start", "prom_end", "down_start", "down_end"]
        ].max(axis=1)

        # check if pos is in any of the metadata entries
        # if so, add all data to results
        for i, row_a in group_chr[1].iterrows():
            for j, row_b in chr_entries.iterrows():
                if row_b["min"] <= row_a["pos"] <= row_b["max"]:
                    # print(row_a)
                    # print(row_b)
                    merged_dict = {**row_a.to_dict(), **row_b.to_dict()}
                    # print(merged_dict)
                    cpg_loc = merged_dict["pos"]
                    merged_dict["CpG location"] = cpg_loc
                    keys = ["chr", "pos", "min", "max"]
                    list(map(merged_dict.pop, keys))
                    merged_dict["CpG region"] = ""
                    if (
                        merged_dict["start pos"] <= cpg_loc <= merged_dict["end pos"]
                    ) or (
                        merged_dict["start pos"] >= cpg_loc >= merged_dict["end pos"]
                    ):
                        merged_dict["CpG region"] = "intragenic"
                    elif (
                        merged_dict["prom_start"] <= cpg_loc <= merged_dict["prom_end"]
                    ) or (
                        merged_dict["prom_start"] >= cpg_loc >= merged_dict["prom_end"]
                    ):
                        merged_dict["CpG region"] = "promoter"
                    elif (
                        merged_dict["down_start"] <= cpg_loc <= merged_dict["down_end"]
                    ) or (
                        merged_dict["down_start"] >= cpg_loc >= merged_dict["down_end"]
                    ):
                        merged_dict["CpG region"] = "downstream"
                    results_as_list_dicts.append(merged_dict)

    return results_as_list_dicts


def make_results_df(results_as_list_dicts):
    results_df = pd.DataFrame(
        results_as_list_dicts,
        columns=[
            "ID",
            "Biomarker name",
            "chrom",
            "start pos",
            "end pos",
            "strand",
            "probe",
            "beta",
            "CpG location",
            "CpG region",
            "N",
            "X",
            "prom_start",
            "prom_end",
            "down_start",
            "down_end",
        ],
    )

    # Change column type to int64 for columns: 'start pos', 'end pos' and 5 other columns
    results_df = results_df.astype(
        {
            "start pos": "int64",
            "end pos": "int64",
            "CpG location": "int64",
            "prom_start": "int64",
            "prom_end": "int64",
            "down_start": "int64",
            "down_end": "int64",
        }
    )

    # Sort by columns: 'id_num' (ascending), 'CpG location' (ascending)
    results_df["id_num"] = results_df["ID"].astype("int64")
    results_df = results_df.sort_values(["id_num", "CpG location"])
    results_df.drop("id_num", axis=1, inplace=True)

    results_df.reset_index(inplace=True, drop=True)

    return results_df


def get_results_df(results_df, panel_data):
    """these steps to merge dfs are broken up because of the long time it takes to go through the dss results"""
    results_as_list_dicts = merge_data_to_list_dicts(
        results_df, panel_data
    )
    results_df = make_results_df(results_as_list_dicts)
    return results_df


def get_region_methylation(data):
    
    # len(data) should > 3
    if len(data) < 3:  # TODO: Discuss the minimum number of modified entries to calculate a region methylated score. Probs will end up using a professional calculation or something from minfi maybe
        raise ValueError("There should be > 3 entries to calculate a region modified score")
    # Find meth/nonmeth for each position in region
    total = 0
    meth = 0
    for cpg_pos in data.iterrows():
        if cpg_pos[1]['beta'] < 0.8:
            total += 1
        else:
            meth += 1
            total += 1
    
    # get region score        
    result = meth/total

    return result, meth, total


def format_results_for_preclin_output(results_df):
    """Make preclin panel output"""
    
    BM_TYPE_FULL = get_BM_TYPE_FULL(path2panel=args.panel)
    
    bm_classif_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)
    panel_result_header_rawmod = preclin_stage_panel_result_header.copy()
    panel_result_header_rawmod.extend(['Meth (beta >= 0.8)', 'Total'])
    bm_classif_panel_rawmod_df = pd.DataFrame(columns=panel_result_header_rawmod)

    for i, (panel_id, data) in enumerate(results_df.groupby(by="ID", observed=False)):
        
        all_mod_data_idxd = all_mod_data.set_index('ID')
        panel_entry = all_mod_data_idxd.loc[panel_id].T.to_dict()
        
        if panel_entry['DNA methylation region'] == 'position':
            # position should only be one entry
            if len(data) > 1:
                raise ValueError("mod position entry has more than one result, this shouldn't happen so need to investigate.")
            print("Check that this line gets the right data, haven't had any testing results available yet.")  #TODO: this
            result = data['beta'][0]
            total = np.nan
            meth = np.nan
            biomarker_type = panel_entry[BM_TYPE_FULL]
            
        elif panel_entry['DNA methylation region'] == 'promoter':
            # calculate promotor region score
            result, meth, total = get_region_methylation(data)
            biomarker_type = panel_entry[BM_TYPE_FULL] + ' - promoter region'
            
        elif panel_entry[BM_TYPE_FULL] == 'expression':
            # calculate promotor region score
            result, meth, total = get_region_methylation(data)
            biomarker_type = panel_entry[BM_TYPE_FULL]
                    
        elif panel_entry[BM_TYPE_FULL] == 'exp_ratio':
            # calculate promotor region score
            result, meth, total = get_region_methylation(data)
            biomarker_type = panel_entry[BM_TYPE_FULL]
            
        else:
            raise ValueError(f"code not ready for DNA methylation region = {panel_entry['DNA methylation region']} or variant type = {panel_entry[BM_TYPE_FULL]}")
        
        if panel_entry[BM_TYPE_FULL] != 'exp_ratio':
            bm_classif_panel_df.loc[i] = [panel_id, panel_entry[BIOMARKER_NAME], panel_entry[SCORING_TYPE], biomarker_type, panel_entry[RESULT_OPTIONS], result] 
        
        bm_classif_panel_rawmod_df.loc[i] = [panel_id, panel_entry[BIOMARKER_NAME], panel_entry[SCORING_TYPE], biomarker_type, panel_entry[RESULT_OPTIONS], result, meth, total] 
        
    return bm_classif_panel_df, bm_classif_panel_rawmod_df


def add_exp_ratio_to_results(bm_classif_panel_df, bm_classif_panel_rawmod_df):
    exp_ratio_panel_results = bm_classif_panel_rawmod_df[bm_classif_panel_rawmod_df['Biomarker Type'] == "exp_ratio"] 
    panel_input_exp_ratio_idxd = panel_data_exp_ratio.set_index('ID')

    exp_ratio_data = dict()
    for panel_id, data in exp_ratio_panel_results.groupby(by="ID", observed=False):
        # print(i)
        data = data.set_index('ID')
        data_entry = data.loc[panel_id].to_dict()
        
        ratio_name = panel_input_exp_ratio_idxd.loc[panel_id]['Expression Ratio Components']
        gene_name = panel_input_exp_ratio_idxd.loc[panel_id][BIOMARKER_NAME]
        region_result = data_entry['Result']
        if ratio_name not in exp_ratio_data.keys():
            exp_ratio_data[ratio_name] = {
                gene_name: (region_result, panel_id)
            }
        else:
            exp_ratio_data[ratio_name][gene_name] = (region_result, panel_id)

    i = len(bm_classif_panel_df)+2
    for k, val in exp_ratio_data.items():
        ratio1, ratio2 = k.split('/')
        
        if ratio1 not in val.keys():
            Warning(f"Expression ratio component {ratio1} not found in exp_ratio_data for {k}, expression ratio unable to be calculated. Consider lowering read depth threshold or sequencing more data.")
            return bm_classif_panel_df
        if ratio2 not in val.keys():
            Warning(f"Expression ratio component {ratio2} not found in exp_ratio_data for {k}, expression ratio unable to be calculated. Consider lowering read depth threshold or sequencing more data.")
            return bm_classif_panel_df
        
        val1, id1 = val[ratio1]
        val2, id2 = val[ratio2]
        print(ratio1, ratio2, val1, val2, id1, id2)
        result = val1/val2
        
        bm_classif_panel_df.loc[i] = [f'{id1}/{id2}', f'{ratio1}/{ratio2}', 'continuous', 'exp_ratio', '0.0-10.0', result]
        i += 1
        
    return bm_classif_panel_df


# Load inputs
with open(args.mod_data, 'r') as fp:
    dss_df = pd.read_csv(fp, sep=',', dtype={'probe': str, 'strand': str})

# Export for methatlas
export_for_methatlas(dss_df, args.out_meth)

# Load panel metadata
panel_data_mod = variant_prep(args.panel, 'mod')
panel_data_exp = variant_prep(args.panel, 'expression')
panel_data_exp_ratio = variant_prep(args.panel, 'exp_ratio')

all_mod_data = pd.concat([panel_data_mod, panel_data_exp, panel_data_exp_ratio])

# Add flanking regions
panel_meth_flank_df = add_prom_start_end(all_mod_data.copy())
panel_meth_flank_df = add_downstream_start_end(panel_meth_flank_df.copy())

# ------------------------------------------------ #
# Compare methylation discovered at those coords to what is expected/useful (from metadata file)
# meth_threshold = 0.8
# `pos + 1` should map correctly to the methylated loci

panel_meth_flank_sorted_df, dss_df_sorted = (
    sort_and_prep_dfs(
        panel_meth_flank_df.copy(), dss_df.copy()
    )
)

# ------------------------------------------------ #
# Create a new dataframe containing all the cpgs in each desired location, getting all info from variants into each entry in results
results_df = get_results_df(
    dss_df_sorted, panel_meth_flank_sorted_df
)

# ------------------------------------------------ #
# Format for panel output
# Additional table that has number of methylated positions in region calculations if needed

bm_classif_panel_df, bm_classif_panel_rawmod_df = format_results_for_preclin_output(results_df)

# Add exp_ratio results to bm_classif_panel_df
bm_classif_panel_df = add_exp_ratio_to_results(bm_classif_panel_df, bm_classif_panel_rawmod_df)

# Minimal extraction of results: align dss positions to panel entries
# This was complex in original; here provide an output structure with placeholder logic
# results_df = pd.DataFrame(columns=["ID","Biomarker name","chrom","start pos","end pos","strand","probe","beta","CpG location","CpG region","N","X","prom_start","prom_end","down_start","down_end"])

# Create placeholder outputs
bm_classif_panel_rawmod_df.to_csv(args.out_raw, index=False)
bm_classif_panel_df.to_csv(args.out_mod, index=False)
