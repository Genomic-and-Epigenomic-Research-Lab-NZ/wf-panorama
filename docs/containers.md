# Containers required for wf-panorama

This document lists the container placeholders used in nextflow.config and module process directives. Build Apptainer images with the described packages or provide equivalent images.

1) general (library/general:latest)
- Python 3, pandas, numpy, cyvcf2, jinja2
- samtools, bedtools, bgzip, tabix, gawk

2) r-methylcibersort (library/r-methylcibersort:latest)
- R 4.x, Bioconductor, minfi, MethylCIBERSORT, ggplot2, other R deps

3) epi2me-wf (library/epi2me-wf:latest)
- nextflow, and dependencies required to run epi2me-labs/wf-human-variation (or rely on the nested pipeline's containers)

4) cibersortx (docker://cibersortx/fractions or local .sif)
- CIBERSORTx runtime container (used via Singularity exec)

Notes:
- For reproducibility, build Apptainer images from Dockerfiles or convert conda envs to container images.
- Once images are built, update nextflow.config params and process.container directives accordingly.
