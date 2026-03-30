nextflow.enable.dsl = 2

workflow samples {
    take: panel_metadata_ch
    emit: samples_ch

    // For now, create an empty samples channel to be populated by template-based sample parsing later.
    Channel.empty().set { samples_ch }
}
