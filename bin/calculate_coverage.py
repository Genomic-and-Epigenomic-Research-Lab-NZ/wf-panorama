#!/usr/bin/env python3
import argparse
import pandas as pd
import os

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


def log_print(message, log):
    """Write message to both the log file and terminal (stdout)."""
    print(message, end='')
    log.write(message)


def calc_coverage(variants_bed, log):
    chr_coverage = dict()
    total_length = 0

    for rowid, entry in variants_bed.iterrows():
        chr_name = entry['#chrom']
        chrStart = entry['chromStart']
        chrEnd = entry['chromEnd']
        if chr_name not in chr_coverage:
            chr_coverage[chr_name] = 0
        chr_coverage[chr_name] += int(chrEnd) - int(chrStart)
        total_length += int(chrEnd) - int(chrStart)

    for chr in CHROMOSOMES:
        if chr not in chr_coverage:
            chr_coverage[chr] = 0
        perc_cov_chr = round((chr_coverage[chr]/CHROM_LENGTHS[chr])*100, 2)
        log_print(f'Length in {chr} total bp (% of chromosome): {chr_coverage[chr]} ({perc_cov_chr}%)\n', log)
        # log_print(f'percent coverage: {(chr_coverage[chr]/CHROM_LENGTHS[chr])*100}\n\n', log)

    perc_HG = round((total_length/HG_LENGTH)*100, 2)
    log_print(f'Number of bp covered in total (% of human genome) {total_length:,} bp ({perc_HG}%)\n', log)
    # log_print(f'percent of HG: {perc_HG}%\n', log)
    return perc_HG


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-bed', required=True)
    parser.add_argument('--output-bed', required=True)
    parser.add_argument('--min-cov', type=float, required=True)
    parser.add_argument('--max-cov', type=float, required=True)
    parser.add_argument('--buffer', type=int, required=True)
    parser.add_argument('--log', required=False, default='coverage.log')
    args = parser.parse_args()

    log = open(args.log, 'w')

    bed_df = pd.read_csv(args.input_bed, sep='\t', names=['#chrom', 'chromStart', 'chromEnd', 'name'])
    perc_HG = calc_coverage(bed_df, log)
    
    log_print(f"Criteria:\n", log)
    log_print(f"  - Buffersize in bp: {args.buffer}\n", log)
    log_print(f"  - Min coverage: {args.min_cov}\n", log)
    log_print(f"  - Max coverage: {args.max_cov}\n", log)

    if args.min_cov < perc_HG < args.max_cov:
        log_print(f"Yay we found the sweet spot! Final coverage: {perc_HG}%\n", log)
        os.system(f'cp {args.input_bed} {args.output_bed}')
    else:
        log_print('Criteria not met. Please adjust buffersize_bp parameter and/or desired min_genome_coverage/max_genome_coverage and rerun. You may use the flag "-resume".\n', log)
        with open(args.output_bed, 'w') as fw:
            fw.write('Criteria not met.')

    log.close()
