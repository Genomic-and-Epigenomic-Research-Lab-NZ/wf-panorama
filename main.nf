#!/usr/bin/env nextflow

// Panorama workflow – DSL2
// Runs sample_processing and then either:
//   - collates results for a BM classifier input (default), or
//   - generates a per-sample clinical report (--clinical_mode true)

nextflow.enable.dsl = 2

include { panel_prep }              from './modules/panel_prep.nf'
include { sample_processing }       from './modules/sample_processing.nf'
include { make_classifier_input }   from './modules/make_classifier_input.nf'
include { get_scores }              from './modules/make_sample_report.nf'
include { generate_report }         from './modules/make_sample_report.nf'

// Aliased publish processes — DSL2 requires a unique process instance per workflow invocation
process publish_target_bed {
    label "wftemplate"
    publishDir ( params.project_outdir, mode: "copy", saveAs: { dirname ? "$dirname/$fname" : fname } )
    input:  tuple path(fname), val(dirname)
    output: path fname
    script: "echo 'Writing output files'"
}

process publish_classifier {
    label "wftemplate"
    publishDir ( params.sample_outdir, mode: "copy" )
    input:  path fname
    output: path fname
    script: "echo 'Writing output files'"
}

process publish_report {
    label "wftemplate"
    publishDir ( params.sample_outdir, mode: "copy" )
    input:  path fname
    output: path fname
    script: "echo 'Writing output files'"
}


// ---------------------------------------------------------------------------
// Named sub-workflows that compose sample_processing with downstream steps
// ---------------------------------------------------------------------------

workflow generate_target_bed {
    take:
        panel_metadata_ch
        publish_dir_ch
    main:
        // Run panel prep to generate target bed for MinKNOW adaptive sampling
        panel_prep(panel_metadata_ch)

        // Publish target bed and fasta for MinKNOW adaptive sampling
        panel_prep.out.minknow_bed
            | map { f -> tuple(f, "minknow_input") } 
            | publish_target_bed
        // panel_prep.out.minknow_fasta
        //     | map { f -> tuple(f, "${params.project_name}/minknow_input") } 
        //     | publish

    emit:
        minknow_bed = panel_prep.out.minknow_bed
        // minknow_fasta = panel_prep.out.minknow_fasta
}

workflow get_classifier_input_sample_data {
    take:
        panel_metadata_ch
        publish_dir_ch
    main:
        // Run per-sample processing (variant calling, methylation, immune, etc.)
        sample_processing(panel_metadata_ch)

        // Collate results for the BM classifier
        make_classifier_input(
            panel_metadata_ch,
            sample_processing.out.snv_panel,
            sample_processing.out.sv_panel,
            sample_processing.out.mod_results,
            sample_processing.out.immune
        )

        // Publish final panel results into params.sample_outdir
        make_classifier_input.out.panel_results
            | publish_classifier

    emit:
        panel_results = make_classifier_input.out.panel_results
}

workflow get_patient_report {
    take:
        panel_metadata_ch
        publish_dir_ch
    main:
        // Run per-sample processing (variant calling, methylation, immune, etc.)
        sample_processing(panel_metadata_ch)

        // Score and generate the clinical report
        get_scores(
            panel_metadata_ch,
            sample_processing.out.snv_panel,
            sample_processing.out.mod_results,
            sample_processing.out.immune
        )
        generate_report(
            panel_metadata_ch,
            Channel.fromPath("${projectDir}/resources/template.md"),
            get_scores.out.scores
        )

        // Publish final report into params.sample_outdir
        generate_report.out.report_md
            .mix(generate_report.out.report_html)
            | publish_report

    emit:
        report = generate_report.out.report_html
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

workflow {
    // Validate required params (all modes)
    if (!params.project_name) {
        error "Please provide --project_name <project_title>"
    }
    if (!params.panel_metadata) {
        error "Please provide --panel_metadata <path_to_csv>"
    }

    // Validate that exactly one mode is specified
    def modeCount = [params.make_target_bed, params.clin_trial_mode, params.clinical_mode].count { it }
    if (modeCount == 0) {
        error "Please specify a mode of operation: --make_target_bed, --clin_trial_mode, or --clinical_mode."
    }
    if (modeCount > 1) {
        error "Please specify only one mode of operation: --make_target_bed, --clin_trial_mode, or --clinical_mode."
    }

    // Validate mode-specific required params
    // TODO: gather all errors into one output to avoid multiple rounds of fixing params and re-running
    def validateSampleProcessingParams = { mode ->
        if (!params.sample)              error "Please provide --sample <sample_name> when using --${mode}"
        if (!params.bam_directory)       error "Please provide --bam_directory <path_to_bam_folder> when using --${mode}"
        if (!params.cibersortx_username) error "Please provide --cibersortx_username <username> when using --${mode}"
        if (!params.cibersortx_token)    error "Please provide --cibersortx_token <token> when using --${mode}"
        if (!params._base_working_dir_user_set) error "Please provide --base_working_dir <path> when using --${mode}. This must point to the same base directory used during --make_target_bed so that shared project files (target bed, panel files) can be located."
    }

    if (params.make_target_bed) {
        // No additional required params for target bed generation
    } else if (params.clin_trial_mode) {
        validateSampleProcessingParams('clin_trial_mode')
    } else if (params.clinical_mode) {
        validateSampleProcessingParams('clinical_mode')
    }

    // Set output directory for results
    def mode_outdir = params.make_target_bed ? params.project_outdir : params.sample_outdir
    publish_dir = file("${mode_outdir}")
    publish_dir.mkdirs()

    // Create input channels
    panel_metadata_ch = Channel.fromPath(params.panel_metadata, checkIfExists: true)
    publish_dir_ch = Channel.fromPath(mode_outdir)

    if (params.make_target_bed) {
        // Generate target bed file for MinKNOW adaptive sampling
        // This process is separate from the main sample processing workflow as it is a different use case
        generate_target_bed(panel_metadata_ch, publish_dir_ch)
    } else if (params.clin_trial_mode) {
        // Generate classifier input for the sample
        get_classifier_input_sample_data(panel_metadata_ch, publish_dir_ch)
    } else if (params.clinical_mode) {
        // Generate individual sample clinical report
        get_patient_report(panel_metadata_ch, publish_dir_ch)
    } else {
        error "Please specify a mode of operation: --make_target_bed, --clin_trial_mode, or --clinical_mode."
    }
}
