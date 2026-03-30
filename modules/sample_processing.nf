nextflow.enable.dsl = 2

include { WF_HUMVAR } from './wf_humvar/wf_humvar.nf'

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
        WF_HUMVAR.out

    // Wire the output from wf_humvar to downstream channels
    wf_humvar_out_ch = WF_HUMVAR.out
    // wf_humvar_out_ch.view { println "wf_humvar output: ${it}" }
}

workflow sample_processing {
    take: panel_metadata_ch
    emit: sample_done_ch

    // Expect params.sample and params.project_name
    def SAMPLE = params.sample
    def PROJECT = params.project_name

    if (!SAMPLE) {
        error 'params.sample must be set to run sample_processing'
    }

    // Input paths
    def bam_pass_dir = params.bam_directory ?: "results/${PROJECT}/${SAMPLE}/bam_pass"
    def reference = file('resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna')
    def targets_bed = file("results/${PROJECT}/minknow_input/targets.bed")
    def tandem_repeat_bed = file('resources/hg38.trf.bed.gz')

    Channel.of(file(bam_pass_dir)).set { bam_dir_ch }
    Channel.of(reference).set { ref_ch }
    Channel.of(targets_bed).set { targets_ch }
    Channel.of(tandem_repeat_bed).set { tr_ch }
    Channel.of(SAMPLE).set { sample_name_ch }
    Channel.of(PROJECT).set { project_name_ch }

    // Call wf_humvar subworkflow
    wf_humvar(bam_dir_ch, ref_ch, targets_ch, tr_ch, sample_name_ch, project_name_ch)

    // Wire the output from wf_humvar to downstream channels
    wf_humvar_out_ch.view { println "wf_humvar output: ${it}" }

    // Map wf_humvar outputs to expected file channels (fallback to known paths)
    def mods1_ch = wf_humvar_out_ch.map { it -> file(it?.mods1 ?: "results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf_mods.1.bedmethyl.gz") }
    def mods2_ch = wf_humvar_out_ch.map { it -> file(it?.mods2 ?: "results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf_mods.2.bedmethyl.gz") }
    def mods_ungrouped_ch = wf_humvar_out_ch.map { it -> file(it?.mods_ungrouped ?: "results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf_mods.ungrouped.bedmethyl.gz") }
    def vcf_clin_ch = wf_humvar_out_ch.map { it -> file(it?.vcf_clinvar ?: "results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf_snp_clinvar.vcf.gz") }
    def vcf_clin_tbi_ch = wf_humvar_out_ch.map { it -> file(it?.vcf_clinvar_tbi ?: "results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf_snp_clinvar.vcf.gz.tbi") }
    def vcf_all_ch = wf_humvar_out_ch.map { it -> file(it?.vcf_all ?: "results/${PROJECT}/${SAMPLE}/wf-humvar/${SAMPLE}.wf_snp.vcf.gz") }

    // combine_bedmethyls
    process combine_bedmethyls {
        tag "combine_bedmethyls.${SAMPLE}"
        cpus 2
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/mod_calling", mode: 'copy'
        input:
            path b1 from mods1_ch
            path b2 from mods2_ch
            path b3 from mods_ungrouped_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.wf_mods.all.bedmethyl.bed" into combined_bed_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}/mod_calling
            bgzip -dc ${b1} ${b2} ${b3} | sort -k1,1 -k2,2n > results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.wf_mods.all.bedmethyl.bed
            """
    }

    process convert_bedmethyl_to_DSS {
        tag "convert_bedmethyl_to_DSS.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/mod_calling", mode: 'copy'
        input:
            path bed from combined_bed_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.wf_mods.all.dss_format.tsv" into dss_ch
        script:
            """
            python3 bin/convert_DSS.py --input ${bed} --output results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.wf_mods.all.dss_format.tsv
            """
    }

    process prep_for_getting_betas {
        tag "prep_for_getting_betas.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/mod_calling", mode: 'copy'
        input:
            path epic from Channel.fromPath('resources/IlluminaEPIC_genomic_locations_hg38.csv')
            path dss from dss_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.pre_beta.csv" into pre_beta_ch
        script:
            """
            python3 bin/process_result_before_betas.py --epic ${epic} --dss ${dss} --out results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.pre_beta.csv
            """
    }

    process add_betas {
        tag "add_betas.${SAMPLE}"
        cpus 1
        memory '16 GB'
        time '2h'
        container params.r_methyl_container ?: "file://${projectDir}/containers/methylcibersort.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/mod_calling", mode: 'copy'
        input:
            path pre from pre_beta_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.post_beta.csv" into post_beta_ch
        script:
            """
            Rscript bin/add_betas.R results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.pre_beta.csv results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.post_beta.csv
            """
    }

    process modification_calling {
        tag "modification_calling.${SAMPLE}"
        cpus 2
        memory '16 GB'
        time '4h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/mod_calling", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path post from post_beta_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.methatlas.csv" into methatlas_ch
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.mod_results.csv" into mod_results_ch
            path "results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.rawmod_results.csv" into rawmod_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}/mod_calling
            python3 bin/modification_calling.py --panel ${panel_meta} --post ${post} --out-meth results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.methatlas.csv --out-mod results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.mod_results.csv --out-raw results/${PROJECT}/${SAMPLE}/mod_calling/${SAMPLE}.rawmod_results.csv
            """
    }

    process run_methylCS {
        tag "run_methylCS.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '1h'
        // container params.r_methyl_container ?: "file://${projectDir}/containers/methylcibersort.sandbox"
        container './containers/methylcibersort.sif'
        // container './containers/methylcibersort.sandbox'
        publishDir "results/${PROJECT}/${SAMPLE}/methylCS", mode: 'copy'
        input:
            path beta from methatlas_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/methylCS/${SAMPLE}.CS_mix_matrix.txt" into cs_mix_ch
            path "results/${PROJECT}/${SAMPLE}/methylCS/${SAMPLE}.CS_bladder_ref.txt" into cs_ref_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}/methylCS
            Rscript bin/methylcibersort.R ${beta} results/${PROJECT}/${SAMPLE}/methylCS/${SAMPLE}.CS_mix_matrix results/${PROJECT}/${SAMPLE}/methylCS/${SAMPLE}.CS_bladder_ref.txt ${SAMPLE}
            """
    }

    process run_CIBERSORTX {
        tag "run_CIBERSORTX.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '2h'
        container params.cibersortx_container ?: 'docker://cibersortx/fractions'
        publishDir "results/${PROJECT}/${SAMPLE}/methylCS", mode: 'copy'
        input:
            path mix from cs_mix_ch
            path sig from cs_ref_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/methylCS/CIBERSORTx_${SAMPLE}_Results.csv" into cibersortx_ch
        script:
            """
            set -euo pipefail
            mkdir -p results/${PROJECT}/${SAMPLE}/methylCS
            cd results/${PROJECT}/${SAMPLE}/methylCS

            MIX_BASENAME=$(basename ${mix})
            SIG_BASENAME=$(basename ${sig})

            # Pass arguments directly to the container ENTRYPOINT (do not call a binary)
            --username ${params.cibersortx_username} --token ${params.cibersortx_token} \
                --mixture /src/data/${MIX_BASENAME} --sigmatrix /src/data/${SIG_BASENAME} \
                --label ${SAMPLE} --perm 1 --QN FALSE --verbose TRUE --outdir /src/outdir
            """
    }

    process run_CIBERSORTX {
        tag "run_CIBERSORTX.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '2h'
        container 'docker://cibersortx/fractions'
        containerOptions = "-B ${mixture.parent}:/src/data -B ${task.workDir}:/src/outdir"
        publishDir "results/${PROJECT}/${SAMPLE}/methylCS", mode: 'copy'

        input:
        path mixture from cs_mix_ch
        path sigmatrix from cs_ref_ch
        val username params.cibersortx_username
        val token params.cibersortx_token
        val sample_name SAMPLE
        val permutations 1 // TODO update num of perms

        output:
        path "outdir"

        script:
        """
        cibersortx_fractions \
            --username ${username} \
            --token ${token} \
            --mixture ${mixture.getName()} \
            --sigmatrix ${sigmatrix.getName()} \
            --label ${sample_name} \
            --perm ${permutations} \
            --QN FALSE \
            --verbose TRUE
        """
    }

    process snv_annotation {
        tag "snv_annotation.${SAMPLE}"
        cpus 2
        memory '8 GB'
        time '2h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/snv_annotation", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path vcf_clin from vcf_clin_ch
            path vcf_clin_tbi from vcf_clin_tbi_ch
            path vcf_all from vcf_all_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/snv_annotation/${SAMPLE}.raw_snv_results.csv" into snv_raw_ch
            path "results/${PROJECT}/${SAMPLE}/snv_annotation/${SAMPLE}.snv_results.csv" into snv_panel_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}/snv_annotation
            python3 bin/snv_annotation.py --panel ${panel_meta} --vcf_clin ${vcf_clin} --vcf_all ${vcf_all} --out_raw results/${PROJECT}/${SAMPLE}/snv_annotation/${SAMPLE}.raw_snv_results.csv --out_panel results/${PROJECT}/${SAMPLE}/snv_annotation/${SAMPLE}.snv_results.csv
            """
    }

    process sv_annotation {
        tag "sv_annotation.${SAMPLE}"
        cpus 2
        memory '8 GB'
        time '2h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/sv_annotation", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path vcf_sv from vcf_all_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/sv_annotation/${SAMPLE}.raw_sv_results.csv" into sv_raw_ch
            path "results/${PROJECT}/${SAMPLE}/sv_annotation/${SAMPLE}.sv_results.csv" into sv_panel_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}/sv_annotation
            python3 bin/sv_annotation.py --panel ${panel_meta} --vcf_sv ${vcf_sv} --out results/${PROJECT}/${SAMPLE}/sv_annotation/${SAMPLE}.raw_sv_results.csv
            """
    }

    process immune_infiltrate_mCS {
        tag "immune_infiltrate_mCS.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}/immune_infiltrate", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path mcs from cibersortx_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/immune_infiltrate/${SAMPLE}.immune_panel_results.csv" into immune_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}/immune_infiltrate
            python3 bin/get_immune_infiltrate.mCS.py --deconv ${mcs} --out results/${PROJECT}/${SAMPLE}/immune_infiltrate/${SAMPLE}.immune_panel_results.csv --panel ${panel_meta}
            """
    }

    process collate_results_for_BM_classifier {
        tag "collate_results.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path snv_res from snv_panel_ch
            path mod_res from mod_results_ch
            path immune_res from immune_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/panel_results.csv" into panel_results_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}
            python3 bin/collate_results_for_BM_classifier.py --panel ${panel_meta} --snv ${snv_res} --mod ${mod_res} --immune ${immune_res} --out results/${PROJECT}/${SAMPLE}/panel_results.csv
            """
    }

    process get_scores {
        tag "get_scores.${SAMPLE}"
        cpus 1
        memory '8 GB'
        time '1h'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "results/${PROJECT}/${SAMPLE}", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path snv_res from snv_panel_ch
            path mod_res from mod_results_ch
            path immune_res from immune_ch
        output:
            path "results/${PROJECT}/${SAMPLE}/${SAMPLE}.scores.csv" into scores_ch
        script:
            """
            mkdir -p results/${PROJECT}/${SAMPLE}
            python3 bin/get_scores.py --panel ${panel_meta} --snv ${snv_res} --mod ${mod_res} --immune ${immune_res} --out results/${PROJECT}/${SAMPLE}/${SAMPLE}.scores.csv
            """
    }

    process generate_report {
        tag "generate_report.${SAMPLE}"
        cpus 1
        memory '4 GB'
        time '30m'
        container "file://${projectDir}/containers/general.sandbox"
        publishDir "workflow/report", mode: 'copy'
        input:
            path panel_meta from Channel.fromPath(params.panel_metadata)
            path template from Channel.fromPath('resources/template.md')
            path scores from scores_ch
        output:
            path "workflow/report/${SAMPLE}.report.md" into report_ch
        script:
            """
            mkdir -p workflow/report
            python3 bin/generate_report.py --panel ${panel_meta} --template ${template} --scores ${scores} --out workflow/report/${SAMPLE}.report.md --sample ${SAMPLE}
            """
    }

    sample_done_ch = report_ch
}
