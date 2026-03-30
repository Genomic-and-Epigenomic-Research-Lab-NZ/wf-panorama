nextflow.enable.dsl = 2

// panel_prep.nf - DSL2 module for panel preparation (resource and container directives added)

workflow panel_prep {
    take: panel_metadata_ch
    emit: panel_bed_ch, all_targets_ch, minknow_bed_ch, minknow_fasta_ch, sorted_targets_ch, ini_targets_ch, final_bed_ch

    // Parameters (from params in main workflow)
    def ref_fasta = file('resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
    def epic_locs = file('resources/IlluminaEPIC_genomic_locations_hg38.csv')
    def immune_ref = file('resources/ref_atlas_bladder.csv')
    def chrom_sizes = file('resources/hg38_no_alt.chrom_sizes')

    // buffer and coverage params from params with defaults
    def buffer_bp = params.buffersize_bp ?: 2000
    def min_cov = params.min_genome_coverage ?: 0.5
    def max_cov = params.max_genome_coverage ?: 2.0
    def project = params.project_name ?: 'project'

    // Create paths for outputs under results/<project>/...
    def panel_bed = file("results/${project}/minknow_input_supp/biomarker_panel.bed")
    def all_targets = file("results/${project}/minknow_input/targets.bed")
    def minknow_bed = file("results/${project}/minknow_input_supp/targets.minknow.${buffer_bp}.bed")
    def minknow_fasta = file("results/${project}/minknow_input_supp/targets.minknow.${buffer_bp}.fasta")
    def sorted_targets = file("results/${project}/minknow_input_supp/sorted_all_targets.${buffer_bp}.bed")
    def ini_targets = file("results/${project}/minknow_input_supp/ini_all_targets.${buffer_bp}.bed")
    def final_bed = file("results/${project}/minknow_input/targets.buffed.bed")

    // Index reference
    Channel.fromPath(ref_fasta.toString())
        .ifEmpty { error "Reference fasta not found: ${ref_fasta}" }
        .set { ref_ch }

    process index_ref {
        tag "index_ref"
        cpus 1
        memory '4 GB'
        time '30m'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "resources", mode: 'copy'

        input:
            path ref from ref_ch
        output:
            path "${ref}.fai" into ref_fai_ch
        script:
            """
            samtools faidx ${ref}
            """
    }

    // Get chrom sizes from .fai
    process get_chrom_sizes {
        tag "get_chrom_sizes"
        cpus 1
        memory '2 GB'
        time '15m'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "resources", mode: 'copy'

        input:
            path fai from ref_fai_ch
        output:
            path "resources/hg38_no_alt.chrom_sizes" into chrom_sizes_ch
        script:
            """
            mkdir -p resources
            cut -f1,2 ${fai} > resources/hg38_no_alt.chrom_sizes
            """
    }

    // Make panel bed
    process make_panel_bed {
        tag "make_panel_bed"
        cpus 2
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${project}/minknow_input_supp", mode: 'copy'

        input:
            path panel_csv from panel_metadata_ch
            path immune_ref_file from Channel.of(immune_ref)
            path epic_locs_file from Channel.of(epic_locs)
        output:
            path panel_bed into panel_bed_out_ch
            path all_targets into all_targets_out_ch
        script:
            """
            mkdir -p $(dirname ${panel_bed})
            mkdir -p $(dirname ${all_targets})
            python3 bin/make_panel_bed.py \
                --panel-csv ${panel_csv} \
                --immune-reference ${immune_ref_file} \
                --epic-locs ${epic_locs_file} \
                --panel-bed ${panel_bed} \
                --all-targets ${all_targets}
            """
    }

    // Run adaptive ref generation
    process run_make_adaptive_ref {
        tag "run_make_adaptive_ref"
        cpus 2
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${project}/minknow_input_supp", mode: 'copy'

        input:
            path targets_bed from all_targets_out_ch
            path ref from ref_ch
            path fai from ref_fai_ch
            path chroms from chrom_sizes_ch
        output:
            path minknow_bed into minknow_bed_ch
            path minknow_fasta into minknow_fasta_ch
            path sorted_targets into sorted_targets_ch
            path ini_targets into ini_targets_ch
        script:
            """
            mkdir -p $(dirname ${minknow_bed})
            mkdir -p $(dirname ${minknow_fasta})
            bash bin/make_adaptive_ref.sh "${targets_bed}" "${ref}" "${fai}" "${chroms}" ${buffer_bp} "${minknow_bed}" "${minknow_fasta}" "${sorted_targets}" "${ini_targets}" ${task.log}
            """
    }

    // Check coverage and produce final buffered targets
    process check_coverage {
        tag "check_coverage"
        cpus 1
        memory '4 GB'
        time '30m'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${project}/minknow_input", mode: 'copy'

        input:
            path minknow_bed_file from minknow_bed_ch
        output:
            path final_bed into final_bed_ch
        script:
            """
            mkdir -p $(dirname ${final_bed})
            python3 bin/calculate_coverage.py --input-bed ${minknow_bed_file} --output-bed ${final_bed} --min-cov ${min_cov} --max-cov ${max_cov} --buffer ${buffer_bp} --log ${task.log}
            """
    }

    // Wire outputs
    panel_bed_out_ch.view { println "panel_bed: ${it}" }
    all_targets_out_ch.view { println "all_targets: ${it}" }

    panel_bed_ch = panel_bed_out_ch
    all_targets_ch = all_targets_out_ch
}
