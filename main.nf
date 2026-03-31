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

// ---------------------------------------------------------------------------
// Named sub-workflows that compose sample_processing with downstream steps
// ---------------------------------------------------------------------------

workflow get_classifier_input_sample_data {
    take:
        panel_metadata_ch
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
    emit:
        panel_results = make_classifier_input.out.panel_results
}

workflow get_patient_report {
    take:
        panel_metadata_ch
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

    // Create input channel
    panel_metadata_ch = Channel.fromPath(params.panel_metadata, checkIfExists: true)

    if (params.get_sample_report) {
        // Generate individual sample clinical report
        get_patient_report(panel_metadata_ch)
    } else {
        // Generate classifier input for the sample
        get_classifier_input_sample_data(panel_metadata_ch)
    }
}
