#!/bin/bash -e
# CLI wrapper derived from make_adaptive_ref.sh
# Usage: make_adaptive_ref.sh <targets_bed> <ref_fasta> <fai> <chrom_sizes> <buffer_bp> <out_minknow_bed> <out_minknow_fasta> <out_sorted> <out_ini> <log>

set -e

TARGETS=${1}
REF=${2}
FAI=${3}
CHROM_SIZES=${4}
BASES_TO_EXPAND_PER_SIDE=${5:-5000}
SLOPPED_BED=${6}
SUBSETTED_FASTA=${7}
SORTED_OUT=${8}
INI_OUT=${9}
LOG=${10:-/dev/stderr}

mkdir -p $(dirname ${SLOPPED_BED})
mkdir -p $(dirname ${SUBSETTED_FASTA})

# sort input BED by chr then start
sort -k 1V,1 -k 2n,2 ${TARGETS} > ${SORTED_OUT}

# bedtools slop
bedtools slop -l ${BASES_TO_EXPAND_PER_SIDE} -r ${BASES_TO_EXPAND_PER_SIDE} -i ${SORTED_OUT} -g ${CHROM_SIZES} > ${INI_OUT}

# merge
bedtools merge -i ${INI_OUT} -c 4 -o collapse > ${SLOPPED_BED}

TOT_WIDTH=$(gawk 'BEGIN{FS="\t"; OFS="\t";tot=0}{tot=tot+$3-$2}END{print tot}' ${SLOPPED_BED})
if [ ${TOT_WIDTH} -lt 500 ]; then
  echo "ERROR: total reference width in ${SLOPPED_BED} is only ${TOT_WIDTH} bps" >> ${LOG}
  exit 1
fi

# extract fasta
bedtools getfasta -fi ${REF} -bed ${SLOPPED_BED} -fo ${SUBSETTED_FASTA} -name

echo "...done" >> ${LOG}
exit 0
