#!/usr/bin/env python3
import argparse
import pandas as pd
from cyvcf2 import VCF
from shared_functions import variant_prep, get_location_string, get_annotation_dict, get_annotation_info_dict, add_result, variant_dict_columns_to_add, preclin_stage_panel_result_header, BIOMARKER_ID, BIOMARKER_NAME, SCORING_TYPE, get_BM_TYPE_FULL, RESULT_OPTIONS

parser = argparse.ArgumentParser()
parser.add_argument('--panel', required=True)
parser.add_argument('--vcf_clin', required=True)
parser.add_argument('--vcf_all', required=True)
parser.add_argument('--out_raw', required=True)
parser.add_argument('--out_panel', required=True)
args = parser.parse_args()

panel_metadata_fp = args.panel
vcf_clinvar_fp = args.vcf_clin
vcf_all_fp = args.vcf_all
snv_output = args.out_raw
snv_preclin_output = args.out_panel

BIOMARKER_TYPE_FULL = get_BM_TYPE_FULL(path2panel=args.panel)

variants_metadata_df_snps = variant_prep(panel_metadata_fp, 'snv')
vcf_clinvar = VCF(vcf_clinvar_fp)
vcf_all = VCF(vcf_all_fp)
vcfs = [vcf_clinvar, vcf_all]

info_to_add_to_metadata = dict()
num_variants_found_total = 0
num_variants_found_clinvar = 0
rows_found = list()
rows_not_found = list()

for i, vcf in enumerate(vcfs):
    clinvar = (i == 0)
    for j, row in enumerate(variants_metadata_df_snps.iterrows()):
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
        for variant in vcf(loc):
            if clinvar:
                info_to_add_to_metadata[target_ID]['ClinVar'] = 'Yes'
                num_variants_found_clinvar += 1
            else:
                info_to_add_to_metadata[target_ID]['ClinVar'] = 'No'
            info_dict = get_annotation_info_dict(info_field = variant.INFO)
            annotation_dict = get_annotation_dict(info_dict_ann = info_dict.get('ANN',''))
            if 'RS' in info_dict.keys():
                dbSNP_ID = info_dict['RS']
                filterdata_dbSNPid = variants_metadata_df_snps.loc[variants_metadata_df_snps['SNP ID'] == f'rs{dbSNP_ID}']
                if len(filterdata_dbSNPid) == 1:
                    info_to_add_to_metadata[target_ID] = add_result(info_to_add_to_metadata[target_ID], variant, annotation_dict, open(args.out_raw + '.log','w'), clinvar=True)
                    entry_found = True
            if entry_found:
                rows_found.append(target_ID)
                break
            if (chrom == variant.CHROM) and (start == variant.start) and (end == variant.end):
                info_to_add_to_metadata[target_ID] = add_result(info_to_add_to_metadata[target_ID], variant, annotation_dict, open(args.out_raw + '.log','w'), clinvar=clinvar)
                entry_found = True
            if entry_found:
                rows_found.append(target_ID)
                break
        if entry_found:
            num_variants_found_total += 1
            rows_found.append(target_ID)
            continue
        else:
            if not clinvar:
                rows_not_found.append(target_ID)

info_df = pd.DataFrame(info_to_add_to_metadata).T
info_df.index.name = 'ID'
info_df.reset_index(inplace=True)
merged_df = pd.merge(variants_metadata_df_snps, info_df, on='ID', how='left')
merged_df.to_csv(snv_output, index=False)

bm_classif_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)
only_genotypes = merged_df[merged_df['Genotype'] != '']
for i, row in only_genotypes.iterrows():
    bm_classif_panel_df.loc[i] = [row[BIOMARKER_ID], row[BIOMARKER_NAME], row[SCORING_TYPE], row[BIOMARKER_TYPE_FULL], row[RESULT_OPTIONS], row['Genotype']]
bm_classif_panel_df.to_csv(snv_preclin_output, index=False)
