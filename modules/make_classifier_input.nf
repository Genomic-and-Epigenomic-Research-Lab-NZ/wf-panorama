process make_classifier_input {
    tag "collate_results.${params.sample}"
    cpus 1
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.sample_outdir}", mode: 'copy'
    input:
        path panel_meta
        path snv_res
        path sv_res
        path mod_res
        path immune_res
    output:
        path "${params.sample}.panel_results.csv", emit: panel_results
    script:
        """
        python3 ${projectDir}/bin/collate_results_for_BM_classifier.py \
            --panel ${panel_meta} \
            --snv ${snv_res} \
            --sv ${sv_res} \
            --mod ${mod_res} \
            --immune ${immune_res} \
            --out_final ${params.sample}.panel_results.csv
        """
}

