#!/bin/bash -e
# Usage: run_epi2me_humvar.sh <bam_dir> <reference> <targets> <tandem_repeat> <sample> <project> <log>

BAM_DIR=${1}
REFERENCE=${2}
TARGETS=${3}
TANDEM_REPEAT=${4}
SAMPLE=${5}
PROJECT=${6}
LOG=${7:-/dev/stderr}

cd results/${PROJECT}/${SAMPLE}/wf-humvar/

nextflow run epi2me-labs/wf-human-variation \
    --bam "${BAM_DIR}" \
    --ref "${REFERENCE}" \
    --bed "${TARGETS}" \
    --out_dir . \
    --sample_name "${SAMPLE}" \
    --sv \
    --snp \
    --mod \
    --str \
    --phased \
    --override_basecaller_cfg "${BASECALLER:-}" \
    --output_gene_summary \
    --output_xam_fmt bam \
    --modkit_args "--preset traditional" \
    --bam_min_coverage 0 \
    --threads ${NXF_THREADS:-4} \
    -profile ${NXF_PROFILE:-local} > ${LOG} 2>&1 || true

touch results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf-humvar.finished.flag
