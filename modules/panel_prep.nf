nextflow.enable.dsl = 2

// panel_prep.nf - DSL2 module for panel preparation

// ---------------------------------------------------------------------------
// Processes (DSL2: defined outside the workflow block; no 'from' / 'into')
// ---------------------------------------------------------------------------

// Note that processes "index_ref" and "get_chrom_sizes" are no longer needed as these files are provided as resources, but they are retained here for completeness and in case they are needed for future reference or modifications
// I want to use them so the user doesn't need to do more stuff, but I am currently not sure how to generate them so they go into the wf-panorama/resources directory rather than the projects directory. Also need to add a check to see if the files already exist in resources and skip these steps if they do, to avoid unnecessary processing.


process index_ref {
    tag "index_ref"
    cpus 1
    memory '4 GB'
    time '30m'
    container "file://${projectDir}/containers/general.sif"
    publishDir "resources", mode: 'copy'

    input:
        path ref

    output:
        path "${ref}.fai", emit: fai

    script:
    """
    samtools faidx ${ref}
    """
}

process get_chrom_sizes {
    tag "get_chrom_sizes"
    cpus 1
    memory '2 GB'
    time '15m'
    container "file://${projectDir}/containers/general.sandbox"
    publishDir "resources", mode: 'copy'

    input:
        path fai

    output:
        path "hg38_no_alt.chrom_sizes", emit: chrom_sizes

    script:
    """
    cut -f1,2 ${fai} > hg38_no_alt.chrom_sizes
    """
}


process get_immune_reference {
    tag "get_immune_reference"
    cpus 1
    memory '2 GB'
    time '15m'
    container params.r_methyl_container ?: "file://${projectDir}/containers/methylcibersort.sif"
    publishDir "${params.out_dir}", mode: 'copy'

    input:
        val cancer_type

    output:
        path "${cancer_type}.mCS_ref.csv", emit: immune_ref

    script:
    """
    micromamba run -n mCS Rscript ${projectDir}/bin/get_mCS_ref_csv.R ${cancer_type} ${cancer_type}.mCS_ref.csv
    """
}


process make_panel_bed {
    tag "make_panel_bed"
    cpus 2
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sif"
    // Only publish the all_targets BED (targets_for_align.bed); suppress biomarker_panel.bed
    publishDir "${params.out_dir}/minknow_input", mode: 'copy',
        saveAs: { filename -> filename == "biomarker_panel.bed" ? null : filename }

    input:
        path panel_csv
        path immune_ref_file
        path epic_locs_file
        val  final_align_bed_name

    output:
        path "biomarker_panel.bed",  emit: panel_bed
        path final_align_bed_name,   emit: all_targets

    script:
    """
    python3 ${projectDir}/bin/make_panel_bed.py \
        --panel-csv ${panel_csv} \
        --immune-reference ${immune_ref_file} \
        --epic-locs ${epic_locs_file} \
        --panel-bed biomarker_panel.bed \
        --all-targets ${final_align_bed_name}
    """
}

process run_make_adaptive_ref {
    tag "run_make_adaptive_ref"
    cpus 2
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sif"
    // No publishDir – outputs are intermediate files only

    input:
        path targets_bed
        path ref
        path fai
        path chroms
        val  buffer_bp

    output:
        path "targets.minknow.${buffer_bp}.bed",    emit: minknow_bed
        path "sorted_all_targets.${buffer_bp}.bed", emit: sorted_targets
        path "ini_all_targets.${buffer_bp}.bed",    emit: ini_targets

    script:
    """
    bash ${projectDir}/bin/make_adaptive_ref.sh \
        "${targets_bed}" \
        "${ref}" \
        "${fai}" \
        "${chroms}" \
        ${buffer_bp} \
        "targets.minknow.${buffer_bp}.bed" \
        "sorted_all_targets.${buffer_bp}.bed" \
        "ini_all_targets.${buffer_bp}.bed"
    """
}
    // "targets.minknow.${buffer_bp}.fasta" \


process check_coverage {
    tag "check_coverage"
    cpus 1
    memory '4 GB'
    time '30m'
    debug true
    container "file://${projectDir}/containers/general.sif"
    // Only publish the final buffered targets BED (targets_buffed.bed)
    publishDir "${params.out_dir}/minknow_input", mode: 'copy'

    input:
        path minknow_bed_file
        val  output_bed_name
        val  min_cov
        val  max_cov
        val  buffer_bp

    output:
        path "${output_bed_name}", emit: final_bed

    script:
    """
    python3 ${projectDir}/bin/calculate_coverage.py \
        --input-bed ${minknow_bed_file} \
        --output-bed ${output_bed_name} \
        --min-cov ${min_cov} \
        --max-cov ${max_cov} \
        --buffer ${buffer_bp}
    """
}

// ---------------------------------------------------------------------------
// panel_prep workflow
// ---------------------------------------------------------------------------

workflow panel_prep {
    take:
        panel_metadata_ch

    main:
        def buffer_bp = params.buffersize_bp       ?: 2000
        def min_cov   = params.min_genome_coverage ?: 0.5
        def max_cov   = params.max_genome_coverage ?: 2.0
        def final_align_bed_name = "${params.project_name}.targets_for_align.bed"
        def final_bed_name       = "${params.project_name}.targets_buffed.bed"

        // def ref_fasta  = file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
        // def ref_fai    = file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai")
        // def epic_locs  = file("${projectDir}/resources/IlluminaEPIC_genomic_locations_hg38.csv")

        def ref_fasta = params.reference_genome 
                        ? file("${projectDir}/resources/${params.reference_genome}") 
                        : file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")

        def ref_fai = params.reference_genome_idx 
                        ? file("${projectDir}/resources/${params.reference_genome_idx}") 
                        : file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai")

        def chrom_sizes_path = params.chrom_sizes_file 
                                ? file("${projectDir}/resources/${params.chrom_sizes_file}") 
                                : file("${projectDir}/resources/hg38_no_alt.chrom_sizes")

        def epic_locs = params.Illumina_epic_locs 
                                ? file("${projectDir}/resources/${params.Illumina_epic_locs}") 
                                : file("${projectDir}/resources/IlluminaEPIC_genomic_locations_hg38.csv")

        ref_ch          = Channel.fromPath(ref_fasta.toString(), checkIfExists: true)
        ref_idx_ch      = Channel.fromPath(ref_fai.toString(), checkIfExists: true)
        chrom_sizes_ch  = Channel.fromPath(chrom_sizes_path.toString(), checkIfExists: true)
        epic_locs_ch    = Channel.fromPath(epic_locs.toString(), checkIfExists: true)
        // immune_ref_ch = Channel.of(immune_ref)

        // 1. Index the reference
        // index_ref(ref_ch)

        // 2. Generate chrom sizes from the .fai
        // get_chrom_sizes(index_ref.out.fai)

        // 1. Generate the immune reference using the R container
        get_immune_reference(Channel.of(params.mCS_cancer_type))

        // 2. Build panel BED and all-targets BED from panel metadata
        make_panel_bed(
            panel_metadata_ch, 
            get_immune_reference.out.immune_ref, 
            epic_locs_ch, 
            final_align_bed_name
        )

        // 3. Generate MinKNOW adaptive-sampling reference files
        run_make_adaptive_ref(
            make_panel_bed.out.all_targets,
            ref_ch,
            ref_idx_ch,
            chrom_sizes_ch,
            Channel.of(buffer_bp)
        )

        // 4. Check coverage and produce final buffered targets BED
        check_coverage(
            run_make_adaptive_ref.out.minknow_bed,
            final_bed_name,
            Channel.of(min_cov),
            Channel.of(max_cov),
            Channel.of(buffer_bp)
        )

    emit:
        panel_bed      = make_panel_bed.out.panel_bed
        all_targets    = make_panel_bed.out.all_targets
        minknow_bed    = run_make_adaptive_ref.out.minknow_bed
        sorted_targets = run_make_adaptive_ref.out.sorted_targets
        ini_targets    = run_make_adaptive_ref.out.ini_targets
        final_bed      = check_coverage.out.final_bed
}
