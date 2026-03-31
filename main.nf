#!/usr/bin/env nextflow

// Panorama workflow – DSL2
// Runs sample_processing and then either:
//   - collates results for a BM classifier input (default), or
//   - generates a per-sample clinical report (--get_sample_report true)

nextflow.enable.dsl = 2

include { sample_processing }    from './modules/sample_processing.nf'
include { make_classifier_input } from './modules/make_classifier_input.nf'
include { get_scores }           from './modules/make_sample_report.nf'
include { generate_report }      from './modules/make_sample_report.nf'


// Use publishDir when possible in the process but this is for when is needed output
// different files. E.g.: outputs from ingress processes or inputs provided by the user.
// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process publish {
    // publish inputs to output directory
    label "wftemplate"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    echo "Writing output files"
    """
}

// ---------------------------------------------------------------------------
// Named sub-workflows that compose sample_processing with downstream steps
// ---------------------------------------------------------------------------

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
            sample_processing.out.mod_results,
            sample_processing.out.immune
        )

        // Publish final panel results into params.out_dir/<sample>/
        make_classifier_input.out.panel_results
            | map { f -> tuple(f, "${params.sample}") }
            | publish

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

        // Publish final report into params.out_dir/<sample>/
        generate_report.out.report
            | map { f -> tuple(f, "${params.sample}") }
            | publish

    emit:
        report = generate_report.out.report
}

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

workflow {
    // Validate required params
    if (!params.sample) {
        error "Please provide --sample <sample_name>"
    }
    if (!params.panel_metadata) {
        error "Please provide --panel_metadata <path_to_csv>"
    }

    // Set output directory for results
    publish_dir = file("${params.out_dir}")
    publish_dir.mkdirs()

    // Create input channels
    panel_metadata_ch = Channel.fromPath(params.panel_metadata, checkIfExists: true)
    publish_dir_ch = Channel.fromPath(params.out_dir)

    if (params.get_sample_report) {
        // Generate individual sample clinical report
        get_patient_report(panel_metadata_ch, publish_dir_ch)
    } else {
        // Generate classifier input for the sample
        get_classifier_input_sample_data(panel_metadata_ch, publish_dir_ch)
    }
}
