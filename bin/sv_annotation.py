#!/usr/bin/env python3
import argparse
import pandas as pd
from cyvcf2 import VCF
from shared_functions import (
    variant_prep, get_location_string, get_annotation_dict, get_annotation_info_dict,
    get_snv_by_genomic_location, add_result, variant_dict_columns_to_add,
    preclin_stage_panel_result_header, BIOMARKER_ID, BIOMARKER_NAME, SCORING_TYPE,
    get_BM_TYPE_FULL, get_full_SCORING_TYPE, RESULT_OPTIONS, RESULT, NOTES
)


def main(args):

    panel_metadata_fp = args.panel
    vcf_sv_fp = args.vcf_sv
    sv_output_raw = args.out_raw
    sv_preclin_out = args.out_panel
    BIOMARKER_TYPE_FULL = get_BM_TYPE_FULL(path2panel=panel_metadata_fp)
    SCORING_TYPE_FULL = get_full_SCORING_TYPE(path2panel=panel_metadata_fp)

    variants_metadata_df_svs = variant_prep(panel_metadata_fp, variant_type='sv')
    vcf_sv = VCF(vcf_sv_fp)

    info_to_add_to_metadata = dict()
    num_variants_found_total = 0
    rows_found = list()
    rows_not_found = list()

    log_fp = args.out_raw + '.log'
    with open(log_fp, 'w') as log:
        for _, row in variants_metadata_df_svs.iterrows():
            target_ID = row['ID']
            if target_ID in rows_found:
                continue
            entry_found = False
            log.write('================\n')
            log.write(f"looking for panel target {target_ID} in sv metadata\n")

            if target_ID not in info_to_add_to_metadata:
                info_to_add_to_metadata[target_ID] = {col: '' for col in variant_dict_columns_to_add}

            loc, chrom, start, end = get_location_string(row_data=row)
            loc_list = [chrom, start, end]

            in_vcf = False
            vars_in_loc = []
            for variant in vcf_sv(loc):
                vars_in_loc.append(variant)
                in_vcf = True
                log.write('------\n')
                info_to_add_to_metadata[target_ID]['ClinVar'] = 'N/A'

                if len(variant.genotypes) > 1:
                    log.write("Multiple sample processing not currently supported.\n")
                    raise ValueError("Multiple sample processing not currently supported.")

                info_dict = get_annotation_info_dict(info_field=variant.INFO)
                annotation_dict = get_annotation_dict(info_dict_ann=info_dict.get('ANN', '')) if info_dict.get('ANN') else {}

                # 1) exact genomic location
                if (chrom == variant.CHROM) and (start == variant.start) and (end == variant.end):
                    info_to_add_to_metadata[target_ID], entry_found = get_snv_by_genomic_location(
                        variants_metadata_df_svs, variant, info_dict, annotation_dict,
                        info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list
                    )
                if entry_found:
                    rows_found.append(target_ID)
                    break

                # 2) end-only genomic location fallback
                if (chrom == variant.CHROM) and (end == variant.end):
                    info_to_add_to_metadata[target_ID], entry_found = get_snv_by_genomic_location(
                        variants_metadata_df_svs, variant, info_dict, annotation_dict,
                        info_to_add_to_metadata[target_ID], entry_found, log, location=loc_list, end_only=True
                    )
                if entry_found:
                    rows_found.append(target_ID)
                    break

            if entry_found:
                num_variants_found_total += 1
                rows_found.append(target_ID)
                
                log.write(f"  target matched: {target_ID}\n\n")
            elif in_vcf:
                log.write("  area is in vcf but precise entry not found\n")
                log.write(str(row.to_dict()) + "\n")
                # TODO: Add variant info to results even if exact variant match was not found.
                if len(vars_in_loc) > 1:
                    log.write(f"  multiple variants found in VCF for {target_ID} at {chrom}:{start}-{end}\n")
                    for v in vars_in_loc:
                        log.write(f"    Variant: {v.CHROM}:{v.start}-{v.end}, INFO: {v.INFO}\n")
                    raise ValueError(f"Multiple variants found in VCF for {target_ID} at {chrom}:{start}-{end}. Not currently supported (feel free to add this functionality).")

                info_dict = get_annotation_info_dict(info_field=variant.INFO)
                annotation_dict = get_annotation_dict(info_dict_ann=info_dict.get('ANN', '')) if info_dict.get('ANN') else {}
                info_to_add_to_metadata[target_ID] = add_result(
                    add_data_to_this_dict=info_to_add_to_metadata[target_ID],
                    variant=variant,
                    annotation_dict=annotation_dict,
                    log=log,
                    note_text="A variant in this approximate location was found but is not an exact match to the required biomarker. Available information is recorded in the sv_annotation logfile."
                )
                result_string = f"Allele: {annotation_dict['Allele']} \nAnnotation: {annotation_dict['Annotation']} \nAnnotation Impact: {annotation_dict['Annotation_Impact']} \nHGVS.c: {annotation_dict['HGVS.c']} \nLength: {info_dict['SVLEN']} \nLocation: {chrom}:{variant.POS}-{variant.end}"
                log.write("  " + result_string + "\n")
            else:
                rows_not_found.append(target_ID)
                log.write(f"  target not found: {target_ID}\n\n")

        log.write(f"number of variants found total: {num_variants_found_total}\n")
        log.write(f"number of variants not found: {len(rows_not_found)}\n")
        log.write(f"number of variants searched: {len(variants_metadata_df_svs)}\n")

    info_df = pd.DataFrame(info_to_add_to_metadata).T
    info_df.index.name = 'ID'
    info_df.reset_index(inplace=True)
    merged_df = pd.merge(variants_metadata_df_svs, info_df, on='ID', how='left')
    
    for i, row in merged_df.iterrows():
        if pd.isna(row[RESULT]):
            merged_df.at[i, RESULT] = 'absent'
        if pd.isna(row[NOTES]):
            merged_df.at[i, NOTES] = 'No variant detected in sample.'

    merged_df.to_csv(sv_output_raw, index=False)

    bm_classif_panel_df = pd.DataFrame(columns=preclin_stage_panel_result_header)

    if len(merged_df) != 0:
        for i, row in merged_df.iterrows():
            bm_classif_panel_df.loc[i] = [row[BIOMARKER_ID], row[BIOMARKER_NAME], row[SCORING_TYPE_FULL], row[BIOMARKER_TYPE_FULL], row[RESULT_OPTIONS], row[RESULT], row[NOTES]]

    bm_classif_panel_df.to_csv(sv_preclin_out, index=False)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--panel', required=True)
    parser.add_argument('--vcf_sv', required=True)
    parser.add_argument('--out_raw', required=True)
    parser.add_argument('--out_panel', required=True)
    args = parser.parse_args()
    
    main(args)

