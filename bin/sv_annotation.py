#!/usr/bin/env python3
import argparse
import pandas as pd
from cyvcf2 import VCF
from shared_functions import variant_prep, get_location_string, get_annotation_dict, get_annotation_info_dict, variant_dict_columns_to_add

parser = argparse.ArgumentParser()
parser.add_argument('--panel', required=True)
parser.add_argument('--vcf_sv', required=True)
parser.add_argument('--out', required=True)
args = parser.parse_args()

panel_metadata_fp = args.panel
vcf_sv_fp = args.vcf_sv
sv_output = args.out

variants_metadata_df_svs = variant_prep(panel_metadata_fp, variant_type='sv')
vcf_sv = VCF(vcf_sv_fp)

info_to_add_to_metadata = dict()
num_variants_found_total = 0
rows_found = list()
rows_not_found = list()

for j, row in enumerate(variants_metadata_df_svs.iterrows()):
    target_ID = row[1]['ID']
    if target_ID in rows_found:
        continue
    entry_found = False
    if target_ID not in info_to_add_to_metadata.keys():
        info_to_add_to_metadata[target_ID] = dict()
        for colname in variant_dict_columns_to_add:
            info_to_add_to_metadata[target_ID][colname] = ''
    loc, chrom, start, end = get_location_string(row_data = row[1])
    loc_list = [chrom, start, end]
    in_vcf = False
    for i, variant in enumerate(vcf_sv(loc)):
        in_vcf = True
        info_to_add_to_metadata[target_ID]['ClinVar'] = 'N/A'
        info_dict = get_annotation_info_dict(info_field = variant.INFO)
        annotation_dict = get_annotation_dict(info_dict_ann = info_dict.get('ANN',''))
        if (chrom == variant.CHROM) and (start == variant.start) and (end == variant.end):
            # Simplified add_result behaviour
            info_to_add_to_metadata[target_ID]['Genotype'] = ''
            entry_found = True
        if entry_found:
            rows_found.append(target_ID)
            break
    if entry_found:
        num_variants_found_total += 1
        rows_found.append(target_ID)
        continue
    elif in_vcf:
        rows_not_found.append(target_ID)
    else:
        rows_not_found.append(target_ID)

info_df = pd.DataFrame(info_to_add_to_metadata).T
info_df.index.name = 'ID'
info_df.reset_index(inplace=True)
merged_df = pd.merge(variants_metadata_df_svs, info_df, on='ID', how='left')
merged_df.to_csv(sv_output, index=False)
