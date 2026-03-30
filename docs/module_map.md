# Module mapping: Snakemake -> Nextflow (Phase 0)

Summary
- DSL2 modules per logical Snakemake rule/group.
- panel_metadata.csv remains a CSV input and is exposed as params.panel_metadata.
- All conda envs will be translated into Apptainer/Singularity container placeholders.

Per-rule mapping (short)

1) run_methylCS
- Nextflow module: modules/run_methylCS.nf
- Inputs: results/{project}/{sample}/mod_calling/{sample}.methatlas.csv
- Outputs: {sample}.CS_mix_matrix.txt, {sample}.CS_bladder_ref.txt
- Script: methylcibersort.r -> wf-panorama/bin/
- Container: methylcibersort.sif (placeholder)
- Notes: requires R + Bioconductor packages

2) run_CIBERSORTX
- Nextflow module: modules/run_CIBERSORTX.nf
- Inputs: mix + sig files
- Outputs: CIBERSORTx_{sample}_Results.csv
- Script: wrapper to call CIBERSORTx via Apptainer or remote image
- Container: cibersortx.sif (placeholder)

3) combine_bedmethyls
- Nextflow module: modules/combine_bedmethyls.nf
- Inputs: three .bedmethyl.gz
- Outputs: .wf_mods.all.bedmethyl.bed
- Script: bgzip -dc ... | sort ...
- Container: general.sif

4) convert_bedmethyl_to_DSS
- Nextflow module: modules/convert_bedmethyl_to_DSS.nf
- Script: convert_DSS.py -> wf-panorama/bin/
- Container: general.sif

5) prep_for_getting_betas
- Nextflow module: modules/prep_for_getting_betas.nf
- Script: process_result_before_betas.py -> wf-panorama/bin/
- Container: general.sif

6) add_betas
- Nextflow module: modules/add_betas.nf
- Script: add_betas.R -> wf-panorama/bin/
- Container: methylcibersort.sif

7) modification_calling
- Nextflow module: modules/modification_calling.nf
- Script: modification_calling.py -> wf-panorama/bin/
- Container: general.sif

8) run_wf_humvar
- Nextflow module: modules/run_wf_humvar.nf
- Script: wrapper to call epi2me-labs/wf-human-variation (external Nextflow)
- Container: use epi2me pipeline containers or call via Apptainer

9) bgzip_clinvar_vcf / tabix_clinvar_vcf_gz
- Nextflow module: modules/bgzip_tabix.nf
- Script: bgzip and tabix
- Container: general.sif

10) snv_annotation
- Nextflow module: modules/snv_annotation.nf
- Script: snv_annotation.py -> wf-panorama/bin/
- Container: general.sif (include cyvcf2)

11) sv_annotation
- Nextflow module: modules/sv_annotation.nf
- Script: sv_annotation.py -> wf-panorama/bin/
- Container: general.sif

12) immune_infiltrate (mCS and methatlas)
- Nextflow modules: modules/immune_infiltrate_mCS.nf and modules/immune_infiltrate_methatlas.nf
- Scripts: get_immune_infiltrate.mCS.py and get_immune_infiltrate.methatlas.py -> wf-panorama/bin/
- Container: general.sif

13) collate_results_for_BM_classifier
- Nextflow module: modules/collate_results_for_BM_classifier.nf
- Script: collate_results_for_BM_classifier.py -> wf-panorama/bin/
- Container: general.sif

14) get_scores
- Nextflow module: modules/get_scores.nf
- Script: get_scores.py -> wf-panorama/bin/
- Container: general.sif

15) generate_report
- Nextflow module: modules/generate_report.nf
- Script: generate_report.py -> wf-panorama/bin/
- Container: general.sif (jinja2)

16) panel_prep group
- Nextflow module(s): modules/panel_prep.nf (subprocesses for index_ref, get_chrom_sizes, make_panel_bed, run_make_adaptive_ref, check_coverage)
- Scripts: make_panel_bed.py, make_adaptive_ref.sh, calculate_coverage.py -> wf-panorama/bin/
- Container: general.sif (samtools, bedtools)

17) merge_bams
- Nextflow module: modules/merge_bams.nf
- Script: merge_bams.sh -> wf-panorama/bin/
- Container: general.sif

18) basecalling (optional)
- Nextflow module: modules/basecalling.nf
- Script: basecalling.sh -> wf-panorama/bin/
- Container: rely on epi2me-labs/wf-basecalling pipeline or build a basecaller image

Special cases & notes
- Nested Nextflow runs (wf-human-variation, wf-basecalling): initial approach is to shell out to those pipelines; later integrate as submodules.
- Python and R scripts using Snakemake objects will be converted to CLI wrappers in wf-panorama/bin/.
- All conda envs will be replaced by Apptainer placeholders; actual image builds are deferred.
- panel_metadata.csv will be passed unchanged as a file param to modules that need it.

Next actions (Phase 1 completed):
- Scaffold and config added (nextflow.config, main.nf) and this module_map created in docs/.
- Phase 2 will implement panel_prep and input handling modules first.
