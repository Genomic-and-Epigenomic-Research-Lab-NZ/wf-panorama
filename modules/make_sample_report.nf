process get_scores {
    tag "get_scores.${params.sample}"
    cpus 1
    memory '8 GB'
    time '1h'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.sample_outdir}", mode: 'copy'
    input:
        path panel_meta
        path snv_res
        path mod_res
        path immune_res
    output:
        path "${params.sample}.scores.csv", emit: scores
        path "${params.sample}.scoring.log", emit: scoring_log
    script:
        """
        python3 ${projectDir}/bin/get_scores.py \
            --panel ${panel_meta} \
            --snv ${snv_res} \
            --mod ${mod_res} \
            --immune ${immune_res} \
            --out ${params.sample}.scores.csv \
            --log "${params.sample}.scoring.log"
        """
}

process generate_report {
    tag "generate_report.${params.sample}"
    cpus 1
    memory '4 GB'
    time '30m'
    container "file://${projectDir}/containers/general.sif"
    publishDir "${params.sample_outdir}", mode: 'copy'
    input:
        path panel_meta
        path template
        path scores
    output:
        path "${params.sample}.report.md", emit: report_md
        path "${params.sample}.report.html", emit: report_html
    script:
        """
        python3 ${projectDir}/bin/generate_report.py \
            --panel ${panel_meta} \
            --template ${template} \
            --scores ${scores} \
            --report_title "${params.report_title}" \
            --project_name ${params.project_name} \
            --sample ${params.sample} \
            --out ${params.sample}.report.md \
            --log "${params.sample}.report.log"
        """
}
