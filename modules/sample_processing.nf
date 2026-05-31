include { WF_HUMVAR } from './wf_humvar/wf_humvar.nf'

// ---------------------------------------------------------------------------
// Processes (DSL2: no 'from' / 'into'; inputs and outputs are positional)
// ---------------------------------------------------------------------------

process combine_bedmethyls {
    tag "combine_bedmethyls.${params.sample}"
    cpus 2
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/mod_calling", mode: 'copy'
    input:
        path b1
        path b2
        path b3
    output:
        path "${params.sample}.wf_mods.all.bedmethyl.bed", emit: combined_bed
    script:
        """
        bgzip -dc ${b1} ${b2} ${b3} | sort -k1,1 -k2,2n > ${params.sample}.wf_mods.all.bedmethyl.bed
        """
}

process convert_bedmethyl_to_DSS {
    tag "convert_bedmethyl_to_DSS.${params.sample}"
    cpus 1
    memory { 8.GB * task.attempt }
    time '1h'
    errorStrategy { task.exitStatus in [137, 138, 139, 140, 143] ? 'retry' : 'terminate' }  // Retry on common OOM exit codes, with increasing memory on each retry
    maxRetries 5
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/mod_calling", mode: 'copy'
    input:
        path bed
    output:
        path "${params.sample}.wf_mods.all.dss_format.tsv", emit: dss
    script:
        """
        python3 ${projectDir}/bin/convert_DSS.py --input ${bed} --output ${params.sample}.wf_mods.all.dss_format.tsv
        """
}

process prep_for_getting_betas {
    tag "prep_for_getting_betas.${params.sample}"
    cpus 1
    memory { 16.GB * task.attempt }
    time { '1h' * task.attempt }
    errorStrategy { task.exitStatus in [137, 138, 139, 140, 143] ? 'retry' : 'terminate' }  // Retry on common OOM exit codes, with increasing memory on each retry
    maxRetries 5
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/mod_calling", mode: 'copy'
    input:
        path epic
        path dss
    output:
        path "${params.sample}.pre_beta.csv", emit: pre_beta
    script:
        """
        python3 ${projectDir}/bin/process_result_before_betas.py --epic ${epic} --dss ${dss} --depth ${params.meth_coverage_threshold} --out ${params.sample}.pre_beta.csv
        """
}

process add_betas {
    tag "add_betas.${params.sample}"
    cpus 1
    memory '16 GB'
    time '2h'
    container params.r_methyl_container ?: "file://${projectDir}/containers/methylcibersort.sif"
    publishDir "${params.out_dir}/${params.sample}/mod_calling", mode: 'copy'
    input:
        path pre
    output:
        path "${params.sample}.post_beta.csv", emit: post_beta
    script:
        """
        micromamba run -n mCS Rscript ${projectDir}/bin/add_betas.R ${pre} ${params.sample}.post_beta.csv
        """
}

process modification_calling {
    tag "modification_calling.${params.sample}"
    cpus 2
    memory '16 GB'
    time '4h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/mod_calling", mode: 'copy'
    input:
        path panel_meta
        path post_betas
    output:
        path "${params.sample}.methatlas.csv",      emit: methatlas
        path "${params.sample}.mod_results.csv",     emit: mod_results
        path "${params.sample}.rawmod_results.csv",  emit: rawmod
    script:
        """
        python3 ${projectDir}/bin/modification_calling.py \
            --panel ${panel_meta} --mod_data ${post_betas} \
            --out-meth ${params.sample}.methatlas.csv \
            --out-mod  ${params.sample}.mod_results.csv \
            --out-raw  ${params.sample}.rawmod_results.csv
        """
}

process run_methylCS {
    tag "run_methylCS.${params.sample}"
    cpus 1
    memory '8 GB'
    time '1h'
    container params.r_methyl_container ?: "file://${projectDir}/containers/methylcibersort.sif"
    publishDir "${params.out_dir}/${params.sample}/immune_infiltrate", mode: 'copy'
    
    input:
        path beta
        val cancer_type
        val sample_name
    
    output:
        path "${params.sample}.CS_mix_matrix.txt", emit: cs_mix
        path "${params.sample}.${cancer_type}.mCS_ref.txt", emit: cs_ref
    
    script:
        """
        micromamba run -n mCS Rscript ${projectDir}/bin/methylcibersort.R ${beta} ${sample_name}.CS_mix_matrix ${params.sample}.${cancer_type}.mCS_ref.txt ${sample_name} ${cancer_type}
        """
}

process run_CIBERSORTX {
    tag "run_CIBERSORTX.${params.sample}"
    cpus 1
    memory '32 GB'
    time '2h'
    container "file://${projectDir}/containers/cibersortx_fractions.sif"
    containerOptions "--bind \${PWD}:/src/data --bind \${PWD}:/src/outdir"
    publishDir "${params.out_dir}/${params.sample}/immune_infiltrate", mode: 'copy', saveAs: { f -> file(f).name }
    input:
        path mixture
        path sigmatrix
        val username
        val token
        val sample_name
        val permutations
    output:
        path "CIBERSORTx_${params.sample}_Results.csv", emit: cibersortx_out
    script:
        """
        /src/CIBERSORTxFractions \\
            --username ${username} \\
            --token ${token} \\
            --mixture ${mixture.getName()} \\
            --sigmatrix ${sigmatrix.getName()} \\
            --label ${sample_name} \\
            --perm ${permutations} \\
            --QN FALSE \\
            --verbose TRUE
        """
}

process snv_prep {
    tag "snv_prep.${params.sample}"
    cpus 1
    memory '4 GB'
    time '30m'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/wf-humvar", mode: 'copy'
    input:
        path vcf_clin_raw
        path vcf_gz
    output:
        path "${params.sample}.wf_snp_clinvar.vcf.gz",     emit: vcf_clin_gz
        path "${params.sample}.wf_snp_clinvar.vcf.gz.tbi", emit: vcf_clin_tbi
        path "${params.sample}.wf_snp.vcf.gz.tbi", emit: vcf_tbi
    script:
        // mv ${vcf_clin_raw}.gz ${params.sample}.wf_snp_clinvar.vcf.gz
        """
        bgzip -k ${vcf_clin_raw}
        tabix ${vcf_clin_raw}.gz
        tabix ${vcf_gz}
        """
}

process snv_annotation {
    tag "snv_annotation.${params.sample}"
    cpus 2
    memory '8 GB'
    time '2h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/snv_annotation", mode: 'copy'
    input:
        path panel_meta
        path vcf_clin
        path vcf_clin_tbi
        path vcf_all
        path vcf_all_tbi
    output:
        path "${params.sample}.raw_snv_results.csv", emit: snv_raw
        path "${params.sample}.snv_results.csv",     emit: snv_panel
    script:
        """
        python3 ${projectDir}/bin/snv_annotation.py \
            --panel ${panel_meta} \
            --vcf_clin ${vcf_clin} \
            --vcf_all ${vcf_all} \
            --out_raw ${params.sample}.raw_snv_results.csv \
            --out_panel ${params.sample}.snv_results.csv
        """
}

process sv_annotation {
    tag "sv_annotation.${params.sample}"
    cpus 2
    memory '8 GB'
    time '2h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/sv_annotation", mode: 'copy'
    input:
        path panel_meta
        path vcf_sv
    output:
        path "${params.sample}.raw_sv_results.csv", emit: sv_raw
        // path "${params.sample}.sv_results.csv",     emit: sv_panel  // TODO: Turned off for now until we get some SV data we can investigate
    script:
        """
        python3 ${projectDir}/bin/sv_annotation.py \
            --panel ${panel_meta} \
            --vcf_sv ${vcf_sv} \
            --out ${params.sample}.raw_sv_results.csv
        """
}

process immune_infiltrate_mCS {
    tag "immune_infiltrate_mCS.${params.sample}"
    cpus 1
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.out_dir}/${params.sample}/immune_infiltrate", mode: 'copy'
    input:
        path panel_meta
        path mcs
    output:
        path "${params.sample}.immune_panel_results.csv", emit: immune
    script:
        """
        python3 ${projectDir}/bin/get_immune_infiltrate.mCS.py \
            --deconv ${mcs} \
            --out ${params.sample}.immune_panel_results.csv \
            --panel ${panel_meta}
        """
}

// ---------------------------------------------------------------------------
// Main sample_processing workflow
// ---------------------------------------------------------------------------

workflow sample_processing {
    take:
        panel_metadata_ch

    main:
        def SAMPLE = params.sample

        if (!SAMPLE) {
            error 'params.sample must be set to run sample_processing'
        }

        // Input channels — use Channel.fromPath for files, Channel.of for scalar values
        bam_dir_ch       = Channel.fromPath(params.bam_directory, type: 'dir', checkIfExists: true)
        // ref_ch           = Channel.fromPath("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna", checkIfExists: true)
        targets_ch       = Channel.fromPath(params.target_bedfile, checkIfExists: true)
        // tr_ch            = Channel.fromPath("${projectDir}/resources/hg38.trf.bed", checkIfExists: true)
        // epic_ch          = Channel.fromPath("${projectDir}/resources/IlluminaEPIC_genomic_locations_hg38.csv", checkIfExists: true)
        sample_name_ch   = Channel.of(SAMPLE)
        project_name_ch  = Channel.of(params.project_name)

        def ref_fasta = params.reference_genome 
                        ? file(params.reference_genome) 
                        : file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
        def ref_fai = params.reference_genome_idx 
                        ? file(params.reference_genome_idx) 
                        : file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai")
        
        def epic_locs = params.Illumina_epic_locs 
                                ? file(params.Illumina_epic_locs) 
                                : file("${projectDir}/resources/IlluminaEPIC_genomic_locations_hg38.csv")

        def ref_trf = params.wf_humvar_trf_file 
                        ? file(params.wf_humvar_trf_file) 
                        : file("${projectDir}/resources/hg38.trf.bed")

        ref_ch          = Channel.fromPath(ref_fasta.toString(), checkIfExists: true)
        ref_fai_ch          = Channel.fromPath(ref_fasta.toString(), checkIfExists: true)
        epic_locs_ch    = Channel.fromPath(epic_locs.toString(), checkIfExists: true)
        tr_ch           = Channel.fromPath(ref_trf.toString(), checkIfExists: true)
        // epic_ch is used by prep_for_getting_betas below


        // 1. Call WF_HUMVAR directly (no passthrough wrapper needed)
        WF_HUMVAR(bam_dir_ch, ref_ch, ref_fai_ch, targets_ch, tr_ch, sample_name_ch, project_name_ch)

        // WF_HUMVAR emits a flat list of files from results/<project>/<sample>/wf-humvar/*
        // Collect them into a single list, then filter by filename pattern.
        wf_humvar_files = WF_HUMVAR.out.humvar_files
            .flatMap { dir -> dir.listFiles().toList() }
            .collect()

        mods1_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_mods\.1\.bedmethyl\.gz$/ }
        }
        mods2_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_mods\.2\.bedmethyl\.gz$/ }
        }
        mods_ungrouped_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_mods\.ungrouped\.bedmethyl\.gz$/ }
        }
        vcf_all_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_snp\.vcf\.gz$/ }
        }
        vcf_clin_raw_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_snp_clinvar\.vcf$/ }
        }

        // 2. Methylation pipeline
        combine_bedmethyls(mods1_ch, mods2_ch, mods_ungrouped_ch)
        convert_bedmethyl_to_DSS(combine_bedmethyls.out.combined_bed)
        prep_for_getting_betas(epic_locs_ch, convert_bedmethyl_to_DSS.out.dss)
        add_betas(prep_for_getting_betas.out.pre_beta)
        modification_calling(panel_metadata_ch, add_betas.out.post_beta)
        run_methylCS(modification_calling.out.methatlas, Channel.of(params.mCS_cancer_type), sample_name_ch)

        // 3. CIBERSORTx deconvolution
        // TODO: remove the skipping if statements and just keep the run block and channel creation 
        if (params.skip_cibersortx) {
            // Use a pre-existing results file from the output directory (e.g. when token is expired)
            cibersortx_out_ch = Channel.fromPath(
                "${params.out_dir}/${params.sample}/immune_infiltrate/CIBERSORTx_${params.sample}_Results.csv",
                checkIfExists: true
            )
        } else {
            run_CIBERSORTX(
                run_methylCS.out.cs_mix,
                run_methylCS.out.cs_ref,
                Channel.of(params.cibersortx_username),
                Channel.of(params.cibersortx_token),
                sample_name_ch,
                Channel.of(params.cibersortx_permutations)
            )
            cibersortx_out_ch = run_CIBERSORTX.out.cibersortx_out
        }

        // 4. SNV / SV annotation
        snv_prep(vcf_clin_raw_ch, vcf_all_ch)
        snv_annotation(panel_metadata_ch, snv_prep.out.vcf_clin_gz, snv_prep.out.vcf_clin_tbi, vcf_all_ch, snv_prep.out.vcf_tbi)
        sv_annotation(panel_metadata_ch, vcf_all_ch)

        // 5. Immune infiltrate
        immune_infiltrate_mCS(panel_metadata_ch, cibersortx_out_ch)

    emit:
        snv_panel    = snv_annotation.out.snv_panel
        // sv_panel     = sv_annotation.out.sv_panel  # TODO: Turned off for now until we get some SV data we can investigate
        sv_panel     = sv_annotation.out.sv_raw
        mod_results  = modification_calling.out.mod_results
        immune       = immune_infiltrate_mCS.out.immune
}
