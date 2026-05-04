#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * WF_HUMVAR subworkflow
 * - Runs epi2me-labs/wf-human-variation inside a process
 * - Emits a channel of files produced under `wf-humvar/**`
 */

workflow WF_HUMVAR {
    take:
        bam_dir_ch
        ref_ch
        targets_ch
        tr_ch
        sample_name_ch
        project_name_ch

    main:
        // Invoke the RUN_WF_HUMVAR process and capture its output channel
        RUN_WF_HUMVAR(
            bam_dir_ch,
            ref_ch,
            targets_ch,
            tr_ch,
            sample_name_ch,
            project_name_ch
        )

    emit:
        // Export the process output channel
        humvar_files = RUN_WF_HUMVAR.out.humvar_files
}

process RUN_WF_HUMVAR {
    tag { "wf-human-variation.${sample_name}" }
    cpus 4
    memory '32 GB'
    time '12h'
    cache 'lenient'  // Cache outputs even if the process script changes (since the script is just a wrapper around a stable wf-human-variation run dir)
    storeDir "${params.out_dir}/${params.sample}/.nextflow_cache"

    // ======== WARNING ========
    // storeDir won't re-run if inputs change but outputs already exist
    // This is a deliberate tradeoff with storeDir — it skips the process purely based on output existence, not input hashes. This process is long, therefore this is intentional. If you change --bed, --bam, etc., the cached wf-humvar directory will still be used.
    // =========================

    input:
        path bam_dir
        path ref
        path targets
        path tr
        val sample_name
        val project_name

    output:
        path "wf-humvar", emit: humvar_files

    script:
    // humvar_run_dir is OUTSIDE the task work dir so it is stable across
    // outer pipeline retries — this is what allows -resume to work on the
    // nested wf-human-variation run.
    def humvar_run_dir = "${params.out_dir}/${params.sample}/wf-humvar-run"
    """
    set -euo pipefail

    # Stable directories for the nested run (outside task workDir)
    mkdir -p \
            ${humvar_run_dir}/.nextflow_home \
            ${humvar_run_dir}/work \
            ${humvar_run_dir}/execution

    # Pin NXF_HOME and NXF_WORK to stable locations so -resume finds the
    # nested pipeline's cache even if the outer task work dir changes.
    export NXF_HOME="${humvar_run_dir}/.nextflow_home"
    export NXF_WORK="${humvar_run_dir}/work"
    export NXF_SINGULARITY_CACHEDIR="${projectDir}/containers/singularity"

    # Output dir for this task (relative, inside task workDir — captured by Nextflow)
    mkdir -p wf-humvar

    echo "Running wf-human-variation for sample ${sample_name} with BAM dir ${bam_dir}"

    nextflow run ${moduleDir}/wf-human-variation \
        --bam ${bam_dir} \
        --ref ${ref} \
        --bed ${targets} \
        --sample_name ${sample_name} \
        --project_name ${project_name} \
        --out_dir wf-humvar \
        --sv \
        --snp \
        --mod \
        --str \
        --phased \
        --output_gene_summary \
        --output_xam_fmt bam \
        --modkit_args "--preset traditional" \
        --bam_min_coverage ${params.wf_humvar_bam_min_coverage} \
        --override_basecaller_cfg "${params.ont_basecaller}" \
        -profile ${params.wf_humvar_profile} \
        -process.executor slurm \
        -w \$NXF_WORK \
        -with-report ${humvar_run_dir}/execution/report.html \
        -with-timeline ${humvar_run_dir}/execution/timeline.html \
        -with-trace ${humvar_run_dir}/execution/trace.txt \
        ${params.wf_humvar_resume ? '-resume' : ''}

    # Ensure outputs are visible to the outer Nextflow
    ls -R wf-humvar
    """
}

// process.executor slurm may not work or be appropriate for running on other systems. Use with caution.
