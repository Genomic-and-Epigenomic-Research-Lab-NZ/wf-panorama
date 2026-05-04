include { WF_HUMVAR } from './wf_humvar/wf_humvar.nf'

// ---------------------------------------------------------------------------
// Processes (DSL2: no 'from' / 'into'; inputs and outputs are positional)
// ---------------------------------------------------------------------------

process combine_bedmethyls {
    tag "combine_bedmethyls.${params.sample}"
    cpus 2
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sandbox"
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
    container "file://${projectDir}/containers/general.sandbox"
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
    container "file://${projectDir}/containers/general.sandbox"
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
    container params.r_methyl_container ?: "file://${projectDir}/containers/methylcibersort.sandbox"
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
    container "file://${projectDir}/containers/general.sandbox"
    publishDir "${params.out_dir}/${params.sample}/mod_calling", mode: 'copy'
    input:
        path panel_meta
        path post
    output:
        path "${params.sample}.methatlas.csv",      emit: methatlas
        path "${params.sample}.mod_results.csv",     emit: mod_results
        path "${params.sample}.rawmod_results.csv",  emit: rawmod
    script:
        """
        python3 ${projectDir}/bin/modification_calling.py \
            --panel ${panel_meta} --post ${post} \
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
    container './containers/methylcibersort.sandbox'
    publishDir "${params.out_dir}/${params.sample}/methylCS", mode: 'copy'
    input:
        path beta
    output:
        path "${params.sample}.CS_mix_matrix.txt", emit: cs_mix
        path "${params.sample}.CS_bladder_ref.txt", emit: cs_ref
    script:
        """
        micromamba run -n mCS Rscript ${projectDir}/bin/methylcibersort.R ${beta} ${params.sample}.CS_mix_matrix ${params.sample}.CS_bladder_ref.txt ${params.sample} ${params.methylcibersort_cancer_type}
        """
}

process run_CIBERSORTX {
    tag "run_CIBERSORTX.${params.sample}"
    cpus 1
    memory '8 GB'
    time '2h'
    container 'docker://cibersortx/fractions'
    publishDir "${params.out_dir}/${params.sample}/methylCS", mode: 'copy'
    input:
        path mixture
        path sigmatrix
        val username
        val token
        val sample_name
    output:
        path "outdir/*", emit: cibersortx_out
    script:
        """
        cibersortx_fractions \
            --username ${username} \
            --token ${token} \
            --mixture ${mixture.getName()} \
            --sigmatrix ${sigmatrix.getName()} \
            --label ${sample_name} \
            --perm 1 \
            --QN FALSE \
            --verbose TRUE
        """
}

process snv_annotation {
    tag "snv_annotation.${params.sample}"
    cpus 2
    memory '8 GB'
    time '2h'
    container "file://${projectDir}/containers/general.sandbox"
    publishDir "${params.out_dir}/${params.sample}/snv_annotation", mode: 'copy'
    input:
        path panel_meta
        path vcf_clin
        path vcf_clin_tbi
        path vcf_all
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
    container "file://${projectDir}/containers/general.sandbox"
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
    container "file://${projectDir}/containers/general.sandbox"
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
// Sub-workflows
// ---------------------------------------------------------------------------

workflow wf_humvar {
    take:
        bam_dir_ch
        ref_ch
        targets_ch
        tr_ch
        sample_name_ch
        project_name_ch

    main:
        WF_HUMVAR(
            bam_dir_ch,
            ref_ch,
            targets_ch,
            tr_ch,
            sample_name_ch,
            project_name_ch
        )

    emit:
        humvar_files = WF_HUMVAR.out.humvar_files
}

// ---------------------------------------------------------------------------
// Main sample_processing workflow
// ---------------------------------------------------------------------------

workflow sample_processing {
    take:
        panel_metadata_ch

    main:
        def SAMPLE  = params.sample
        def PROJECT = params.project_name

        if (!SAMPLE) {
            error 'params.sample must be set to run sample_processing'
        }

        // Input paths
        def bam_pass_dir     = file(params.bam_directory)
        def reference        = file("${projectDir}/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna")
        def targets_bed      = file(params.target_bedfile)
        // def tandem_repeat_bed = file("${projectDir}/resources/hg38.trf.bed.gz")
        def tandem_repeat_bed = file("${projectDir}/resources/hg38.trf.bed")
        def epic_file        = file("${projectDir}/resources/IlluminaEPIC_genomic_locations_hg38.csv")

        Channel.of(bam_pass_dir).set       { bam_dir_ch }
        Channel.of(reference).set          { ref_ch }
        Channel.of(targets_bed).set        { targets_ch }
        Channel.of(tandem_repeat_bed).set  { tr_ch }
        Channel.of(SAMPLE).set             { sample_name_ch }
        Channel.of(PROJECT).set            { project_name_ch }

        // 1. Call wf_humvar sub-workflow
        wf_humvar(bam_dir_ch, ref_ch, targets_ch, tr_ch, sample_name_ch, project_name_ch)

        // WF_HUMVAR emits a flat list of files from results/<project>/<sample>/wf-humvar/*
        // Collect them into a single list, then filter by filename pattern.
        wf_humvar_files = wf_humvar.out.humvar_files
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
        vcf_clin_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_snp_clinvar\.vcf\.gz$/ }
        }
        vcf_clin_tbi_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_snp_clinvar\.vcf\.gz\.tbi$/ }
        }
        vcf_all_ch = wf_humvar_files.map { files ->
            files.find { it.name =~ /\.wf_snp\.vcf\.gz$/ }
        }

        // 2. Methylation pipeline
        combine_bedmethyls(mods1_ch, mods2_ch, mods_ungrouped_ch)
        convert_bedmethyl_to_DSS(combine_bedmethyls.out.combined_bed)
        prep_for_getting_betas(Channel.of(epic_file), convert_bedmethyl_to_DSS.out.dss)
        add_betas(prep_for_getting_betas.out.pre_beta)
        modification_calling(panel_metadata_ch, add_betas.out.post_beta)
        run_methylCS(modification_calling.out.methatlas)

        // 3. CIBERSORTx deconvolution
        run_CIBERSORTX(
            run_methylCS.out.cs_mix,
            run_methylCS.out.cs_ref,
            params.cibersortx_username,
            params.cibersortx_token,
            SAMPLE
        )

        // 4. SNV / SV annotation
        snv_annotation(panel_metadata_ch, vcf_clin_ch, vcf_clin_tbi_ch, vcf_all_ch)
        sv_annotation(panel_metadata_ch, vcf_all_ch)

        // 5. Immune infiltrate
        immune_infiltrate_mCS(panel_metadata_ch, run_CIBERSORTX.out.cibersortx_out)

    emit:
        snv_panel    = snv_annotation.out.snv_panel
        // sv_panel     = sv_annotation.out.sv_panel
        sv_panel     = sv_annotation.out.sv_raw
        mod_results  = modification_calling.out.mod_results
        immune       = immune_infiltrate_mCS.out.immune
}
