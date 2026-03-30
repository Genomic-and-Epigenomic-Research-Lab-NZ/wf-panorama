Make the CLI scripts executable:

cd /projects/uow/GERL/dejlu879/Panorama/wf-panorama/bin
chmod +x *.py *.sh *.r

Notes:
- The pipeline uses Apptainer/Singularity container placeholders defined in nextflow.config and module files.
- Build or provide appropriate Singularity images for: general, r-methylcibersort, epi2me-wf, cibersortx, etc.

