#!/bin/bash -e
# CLI wrapper derived from make_adaptive_ref.sh
# ----------------------- #
# Original script by
# Stephane Plaisance (VIB-NC) 2021/03/23; v1.0
# visit our Git: https://github.com/Nucleomics-VIB

# Modified by Lucy Picard for Panorama
# ----------------------- #

# version="1.1, 2025_02_08 - LP"
# version="1.2, 2025_10_03 - LP"
# version="1.3, 2026_04_20 - LP"
version="1.4, 2026_08_03 - LP"
# Usage: make_adaptive_ref.sh <targets_bed> <ref_fasta> <fai> <chrom_sizes> <buffer_bp> <out_minknow_bed> <out_sorted> <out_ini> <log>

set -e

TARGETS=${1}
REF=${2}
FAI=${3}
CHROM_SIZES=${4}
BASES_TO_EXPAND_PER_SIDE=${5:-5000}
SLOPPED_BED=${6}
# SUBSETTED_FASTA=${7}
SORTED_OUT=${7}
INI_OUT=${8}
LOG=${9:-/dev/stderr}

mkdir -p $(dirname ${SLOPPED_BED})
# mkdir -p $(dirname ${SUBSETTED_FASTA})

# sort input BED by chr then start
echo "Sorting input BED file ${TARGETS} by chr then start" >> ${LOG}
sort -k 1V,1 -k 2n,2 ${TARGETS} > ${SORTED_OUT}

# bedtools slop
echo "Expanding ${SORTED_OUT} by ${BASES_TO_EXPAND_PER_SIDE} bps on each side" >> ${LOG}
bedtools slop -l ${BASES_TO_EXPAND_PER_SIDE} -r ${BASES_TO_EXPAND_PER_SIDE} -i ${SORTED_OUT} -g ${CHROM_SIZES} > ${INI_OUT}

# merge
echo "Merging overlapping regions in ${INI_OUT}" >> ${LOG}
bedtools merge -i ${INI_OUT} -c 4 -o collapse > ${SLOPPED_BED}

# check that the total width of the reference is at least 500 bps
echo "Checking that the total width of the reference in ${SLOPPED_BED} is at least 500 bps" >> ${LOG}
TOT_WIDTH=$(gawk 'BEGIN{FS="\t"; OFS="\t";tot=0}{tot=tot+$3-$2}END{print tot}' ${SLOPPED_BED})
if [ ${TOT_WIDTH} -lt 500 ]; then
  echo "ERROR: total reference width in ${SLOPPED_BED} is only ${TOT_WIDTH} bps" >> ${LOG}
  exit 1
fi

# extract fasta
# bedtools getfasta -fi ${REF} -bed ${SLOPPED_BED} -fo ${SUBSETTED_FASTA} -name

echo "...done" >> ${LOG}
exit 0
