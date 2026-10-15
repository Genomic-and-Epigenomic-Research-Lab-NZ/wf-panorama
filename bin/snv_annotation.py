#!/usr/bin/env python3
import argparse
import pandas as pd
from cyvcf2 import VCF
from shared_functions import (
    variant_prep, get_location_string, get_annotation_dict, get_annotation_info_dict,
    add_result, get_snv_by_genomic_location, variant_dict_columns_to_add,
    preclin_stage_panel_result_header, BIOMARKER_ID, BIOMARKER_NAME, SCORING_TYPE,
    get_BM_TYPE_FULL, get_full_SCORING_TYPE, RESULT_OPTIONS, RESULT, NOTES
)


def get_snv_by_RS_number(variants_metadata_df_sub, variant, info_dict, annotation_dict, row_data, entry_found, log, clinvar=False):
    dbSNP_ID = info_dict['RS']
    filterdata_dbSNPid = variants_metadata_df_sub.loc[variants_metadata_df_sub['SNP ID'] == f'rs{dbSNP_ID}']
    if len(filterdata_dbSNPid) == 1:
        log.write("  entry found by dbSNP ID\n")
        row_data = add_result(row_data, variant, annotation_dict, log, clinvar=clinvar)
        entry_found = True
    elif len(filterdata_dbSNPid) > 1:
        log.write(f"  ERROR: ambiguous dbSNP match for rs{dbSNP_ID}\n")
        raise ValueError(f"Ambiguous dbSNP match for rs{dbSNP_ID}")
    return row_data, entry_found

def main(args):

    panel_metadata_fp = args.panel
    vcf_clinvar_fp = args.vcf_clin
    vcf_all_fp = args.vcf_all
    snv_output = args.out_raw
    snv_preclin_output = args.out_panel

    BIOMARKER_TYPE_FULL = get_BM_TYPE_FULL(path2panel=panel_metadata_fp)
    SCORING_TYPE_FULL = get_full_SCORING_TYPE(path2panel=panel_metadata_fp)

    variants_metadata_df_snps = variant_prep(panel_metadata_fp, 'snv')
    vcf_clinvar = VCF(vcf_clinvar_fp)
    vcf_all = VCF(vcf_all_fp)
    vcfs = [vcf_clinvar, vcf_all]

    info_to_add_to_metadata = dict()
    num_variants_found_total = 0
    num_variants_found_clinvar = 0
    rows_found = list()
    rows_not_found = list()

    log_fp = args.out_raw + '.log'
    with open(log_fp, 'w') as log:
        for i, vcf in enumerate(vcfs):
            clinvar = (i == 0)
            log.write("################## CLINVAR ######################\n" if clinvar else "################## NOT CLINVAR ##################\n")

            for _, row in variants_metadata_df_snps.iterrows():
                target_ID = row['ID']
                if target_ID in rows_found:
                    continue

                log.write(f"looking for panel target {target_ID} in snv metadata\n")
                entry_found = False
                if target_ID not in info_to_add_to_metadata:
                    info_to_add_to_metadata[target_ID] = {col: '' for col in variant_dict_columns_to_add}

                loc, chrom, start, end = get_location_string(row_data=row)
                loc_list = [chrom, start, end]

                for variant in vcf(loc):
                    log.write("================\n")
                    if len(variant.genotypes) > 1:
                        log.write("Multiple sample processing not currently supported.\n")
                        raise NotImplementedError("Multiple sample processing not currently supported.")

                    info_dict = get_annotation_info_dict(info_field=variant.INFO)
                    annotation_dict = get_annotation_dict(info_dict_ann=info_dict.get('ANN', '')) if info_dict.get('ANN') else {}

                    # 1) RS match
                    if 'RS' in info_dict:
                        info_to_add_to_metadata[target_ID], entry_found = get_snv_by_RS_number(
                            variants_metadata_df_snps, variant, info_dict, annotation_dict,
                            info_to_add_to_metadata[target_ID], entry_found, log, clinvar=clinvar
                        )
                    if entry_found:
                        info_to_add_to_metadata[target_ID]['ClinVar'] = 'Yes' if clinvar else 'No'
                        if clinvar:
                            num_variants_found_clinvar += 1
                        break

                    # 2) exact genomic match
                    if (chrom == variant.CHROM) and (start == variant.start) and (end == variant.end):
                        info_to_add_to_metadata[target_ID], entry_found = get_snv_by_genomic_location(
                            variants_metadata_df_snps, variant, info_dict, annotation_dict,
                            info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list, clinvar=clinvar
                        )
                    if entry_found:
                        info_to_add_to_metadata[target_ID]['ClinVar'] = 'Yes' if clinvar else 'No'
                        if clinvar:
                            num_variants_found_clinvar += 1
                        break

                    # 3) end-only fallback
                    if (chrom == variant.CHROM) and (end == variant.end):
                        info_to_add_to_metadata[target_ID], entry_found = get_snv_by_genomic_location(
                            variants_metadata_df_snps, variant, info_dict, annotation_dict,
                            info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list, end_only=True, clinvar=clinvar
                        )
                    if entry_found:
                        info_to_add_to_metadata[target_ID]['ClinVar'] = 'Yes' if clinvar else 'No'
                        if clinvar:
                            num_variants_found_clinvar += 1
                        break

                    # 4) start-only fallback - likely indel
                    if (chrom == variant.CHROM) and (start == variant.start):
                        info_to_add_to_metadata[target_ID], entry_found = get_snv_by_genomic_location(
                            variants_metadata_df_snps, variant, info_dict, annotation_dict,
                            info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list, start_only=True, clinvar=clinvar
                        )
                    if entry_found:
                        info_to_add_to_metadata[target_ID]['ClinVar'] = 'Yes' if clinvar else 'No'
                        if clinvar:
                            num_variants_found_clinvar += 1
                        break

                if entry_found:
                    num_variants_found_total += 1
                    rows_found.append(target_ID)
                    log.write(f"  target matched: {target_ID}\n\n")
                elif not clinvar:
                    rows_not_found.append(target_ID)
                    log.write(f"  target not found after non-ClinVar pass: {target_ID}\n\n")

            log.write(f"\nnumber of variants found total: {num_variants_found_total}\n")
            log.write(f"  number of variants found in clinvar: {num_variants_found_clinvar}\n")
            log.write(f"number of variants not found: {len(rows_not_found)}\n")
            log.write(f"number of variants searched: {len(variants_metadata_df_snps)}\n")

    info_df = pd.DataFrame(info_to_add_to_metadata).T
    info_df.index.name = 'ID'
    info_df.reset_index(inplace=True)
    merged_df = pd.merge(variants_metadata_df_snps, info_df, on='ID', how='left')
    merged_df.to_csv(snv_output, index=False)

    bm_classif_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)
    if len(merged_df) != 0:
        only_genotypes = merged_df[merged_df['Genotype'] != '']
        for i, row in only_genotypes.iterrows():
            bm_classif_panel_df.loc[i] = [row[BIOMARKER_ID], row[BIOMARKER_NAME], row[SCORING_TYPE_FULL], row[BIOMARKER_TYPE_FULL], row[RESULT_OPTIONS], row['Genotype'], row[NOTES]]

    bm_classif_panel_df.to_csv(snv_preclin_output, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel', required=True)
    parser.add_argument('--vcf_clin', required=True)
    parser.add_argument('--vcf_all', required=True)
    parser.add_argument('--out_raw', required=True)
    parser.add_argument('--out_panel', required=True)
    args = parser.parse_args()
    
    main(args)
