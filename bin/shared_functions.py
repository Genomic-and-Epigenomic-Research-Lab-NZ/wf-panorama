#!/usr/bin/env python3
# Local copy of shared_functions.py adapted for wf-panorama
import pandas as pd

HG_LENGTH = 3100000000

CHROMOSOMES = [
    "chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
    "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
    "chr21","chr22","chrX","chrY","chrM",
]

CHROM_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
    "chrM": 16569
}


def get_BM_TYPE_FULL(path2panel):
    """
    Ensures we always get the updated list that is provided to the user
    """
    with open(path2panel, "r") as f:
        header = f.readline().strip('\ufeff').strip('\n')
    header_split = header.split('"')
    for h in header_split:
        if h.startswith('Biomarker Type'):  # Biomarker Type header
            return h

def get_VARIANT_TYPES(BM_TYPE_FULL):
    return BM_TYPE_FULL.split('(')[1].split(')')[0].split(', ')

def get_full_SCORING_TYPE(path2panel):
    """
    Ensures we always get the updated list that is provided to the user
    """
    with open(path2panel, "r") as f:
        header = f.readline().strip('\ufeff').strip('\n')
    header_split = header.split('"')
    for h in header_split:
        if h.startswith('Scoring'):  # Scoring Type header
            return h

# BIOMARKER_TYPE_FULL = "Biomarker Type (snv, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, immune_inf, microsatellite, demographic, clinicopathology)"
# VARIANT_TYPE = "Biomarker Type (snv, sv, mod, area_mutations, expression, exp_ratio, immune_ratio, immune_inf, microsatellite, demographic, clinicopathology)" 

# BIOMARKER_TYPE_FULL is no longer set at module level.
# Each script that needs it should call: BIOMARKER_TYPE_FULL = get_BM_TYPE_FULL(path2panel=args.panel)
# VARIANT_TYPE should no longer be used, and BIOMARKER_TYPE_FULL should be used instead

BIOMARKER_ID = "ID"
BIOMARKER_NAME = "Biomarker name"
BIOMARKER_TYPE = "Biomarker Type"
SCORING_TYPE = "Scoring Type"
RESULT_OPTIONS = "Result Options"
RESULT = "Result"
NOTES = "Notes"
EXP_RATIO_COMPONENTS = "Expression Ratio Components"

preclin_stage_panel_result_header = [BIOMARKER_ID, BIOMARKER_NAME, SCORING_TYPE, BIOMARKER_TYPE, RESULT_OPTIONS, RESULT, NOTES]

variant_dict_columns_to_add = ['ClinVar', 'Significance (ClinVar)', 'Consequence (Clinvar)', 'Reference Allele', 'Variant Allele', 'Genotype', 'HGVS.c', 'HGVS.p', 'SV Length', 'SV Type', NOTES]


def filter_to_variant_type(panel_metadata_df, variant_type, BM_TYPE_FULL):
    VARIANT_TYPES = get_VARIANT_TYPES(BM_TYPE_FULL)
    if variant_type not in VARIANT_TYPES:
        raise ValueError("variant_type must be one of: " + ",".join(VARIANT_TYPES))
    panel_metadata_df = panel_metadata_df[(panel_metadata_df[BM_TYPE_FULL] == variant_type)]
    return panel_metadata_df


def reduce_metadata_df(df):
    cols_to_drop = ["Is variant in coding region?","Illumina EPIC ID","Notes","References","length"]
    return df.drop(columns=[c for c in cols_to_drop if c in df.columns], errors='ignore')


def variant_prep(path2_variants_metadata_csv, variant_type):
    with open(path2_variants_metadata_csv, "r") as fpc:
        variants_metadata_df = pd.read_csv(fpc, dtype={"ID": str})
    columns_to_drop = variants_metadata_df.filter(like="Unnamed").columns
    variants_metadata_df.drop(columns=columns_to_drop, inplace=True)
    BM_TYPE_FULL = get_BM_TYPE_FULL(path2panel=path2_variants_metadata_csv)
    variants_metadata_df_v = filter_to_variant_type(variants_metadata_df, variant_type=variant_type, BM_TYPE_FULL=BM_TYPE_FULL)
    variants_metadata_df_v = reduce_metadata_df(variants_metadata_df_v)
    return variants_metadata_df_v


def get_location_string(row_data):
    chrom = row_data['chrom']
    start = int(str(row_data['start pos']).replace(',', ''))
    end = int(str(row_data['end pos']).replace(',', ''))
    loc = (f'{chrom}:{start - 2}-{end + 2}')
    return (loc, chrom, start, end)


def get_annotation_dict(info_dict_ann):
    annotation_header = 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
    annotation_dict = dict(zip(annotation_header.split(' | '), info_dict_ann.split('|')))
    return annotation_dict


def get_annotation_info_dict(info_field):
    info_dict = {}
    for entry in info_field:
        info_dict[entry[0]] = entry[1]
    return info_dict


def add_result(add_data_to_this_dict, variant, annotation_dict, log, clinvar=False, note_text=None):
    from warnings import warn
    if not variant.FILTERS[0] == 'PASS':
        log.write(f"  WARNING: This variant did not pass filters. Variant position: {variant.CHROM}:{variant.start}-{variant.end}\n")
        warn(f"This variant did not pass filters. Variant position: {variant.CHROM}:{variant.start}-{variant.end}")
    ref_allele = variant.REF
    if len(variant.genotypes) > 1:
        log.write("Multiple sample processing not supported.\n")
        raise ValueError("Multiple sample processing not supported.")
    genotype_access_list = [variant.REF]
    genotype_access_list.extend(variant.ALT)
    for genotype in variant.genotypes:
        allele1, allele2, phased = genotype
        allele_base_1 = genotype_access_list[allele1]
        allele_base_2 = genotype_access_list[allele2]
        genotype_text = (f"{allele_base_1}|{allele_base_2}")
    try:
        annotation_consequence = annotation_dict.get('Annotation', 'N/Av')
    except:
        annotation_consequence = 'N/Av'
    try:
        significance = variant.INFO.get('CLNSIG')
    except:
        significance = 'N/Av'
    try:
        sv_length = variant.INFO.get('SVLEN')
    except:
        sv_length = 'N/A'
    try:
        sv_type = variant.INFO.get('SVTYPE')
    except:
        sv_type = 'N/A'
    
    try:
        location = f"{variant.CHROM}:{variant.POS}-{variant.end}"
    except:
        location = 'N/A'
        
    if sv_type != 'N/A':
        try:
            allele = annotation_dict.get('Allele', 'N/A')
        except:
            allele = 'N/A'
        try:
            annotation_impact = annotation_dict.get('Annotation_Impact', 'N/A')
        except:
            annotation_impact = 'N/A'
    else:
        allele = 'N/A'
        annotation_impact = 'N/A'
        
    
    annotation_HGVSc = annotation_dict.get('HGVS.c', '')
    annotation_HGVSp = annotation_dict.get('HGVS.p', '')
    add_data_to_this_dict['Significance (ClinVar)'] = significance
    add_data_to_this_dict['Consequence (Clinvar)'] = annotation_consequence
    add_data_to_this_dict['Reference Allele'] = ref_allele
    add_data_to_this_dict['Genotype'] = genotype_text
    add_data_to_this_dict['HGVS.c'] = annotation_HGVSc
    add_data_to_this_dict['HGVS.p'] = annotation_HGVSp
    add_data_to_this_dict['SV Length'] = sv_length
    add_data_to_this_dict['SV Type'] = sv_type
    add_data_to_this_dict['Allele'] = sv_type
    add_data_to_this_dict['Annotation Impact'] = annotation_impact
    add_data_to_this_dict['Location'] = location
    add_data_to_this_dict['Notes'] = note_text if note_text is not None else ''
    
    return add_data_to_this_dict


def get_snv_by_genomic_location(
    variants_metadata_df_sub,
    variant,
    info_dict,
    annotation_dict,
    row_data,
    entry_found,
    log,
    location=list(),
    end_only=False,
    start_only=False,
    clinvar=False
):
    """
    Match SNV by genomic location using:
    - exact chrom/start/end
    - end-only fallback
    - start-only fallback (indel-like)
    """
    chrom, start, end = location

    chrom_col = variants_metadata_df_sub['chrom'].astype(str)
    start_col = variants_metadata_df_sub['start pos'].astype(str).str.replace(',', '', regex=False)
    end_col = variants_metadata_df_sub['end pos'].astype(str).str.replace(',', '', regex=False)

    if end_only:
        filterdata_loc = variants_metadata_df_sub.loc[
            (chrom_col == str(chrom)) & (end_col == str(end))
        ]
    elif start_only:
        filterdata_loc = variants_metadata_df_sub.loc[
            (chrom_col == str(chrom)) & (start_col == str(start))
        ]
    else:
        filterdata_loc = variants_metadata_df_sub.loc[
            (chrom_col == str(chrom)) & (start_col == str(start)) & (end_col == str(end))
        ]

    if len(filterdata_loc) == 1:
        if end_only:
            log.write("  entry found by end-only genomic location\n")
        elif start_only:
            log.write("  entry found by start-only genomic location (possible indel)\n")
        else:
            log.write("  entry found by exact genomic location\n")
        row_data = add_result(
            add_data_to_this_dict=row_data,
            variant=variant,
            annotation_dict=annotation_dict,
            log=log,
            clinvar=clinvar
        )
        entry_found = True
    elif len(filterdata_loc) > 1:
        log.write(f"  ERROR: ambiguous genomic match at {chrom}:{start}-{end} ({len(filterdata_loc)} rows)\n")
        raise ValueError(f"Ambiguous genomic match at {chrom}:{start}-{end}")

    return row_data, entry_found
