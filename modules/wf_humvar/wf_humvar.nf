nextflow.enable.dsl = 2

process WF_HUMVAR {
    tag { "${sample_name}" }
    cpus 4
    memory '32 GB'
    time '12h'
    container 'epi2me-labs/wf-human-variation:sha2b856c1f358ddf1576217a336bc0e9864b6dc0ed'

    input:
        path bam_dir
        path ref
        path targets
        path tr
        val sample_name
        val project_name

    output:
        path "results/${project_name}/${sample_name}/wf-humvar/*"

    script:
    """
    mkdir -p results/${project_name}/${sample_name}/wf-humvar

    wf-human-variation --bam ${bam_dir} \
            --ref ${ref} \
            --bed ${targets} \
            --tr ${tr} \
            --sample_name ${sample_name} \
            --project_name ${project_name} \
            --out_dir results/${project_name}/${sample_name}/wf-humvar \
            --threads ${params.wf_humvar_threads ?: 8}
    """
}

// Old nested workflow

// nextflow.enable.dsl = 2

// workflow wf_humvar {
//     take: bam_dir_ch, ref_ch, targets_ch, tr_ch, sample_name_ch, project_name_ch
//     emit: wf_humvar_out_ch

//     process exec_humvar {
//         tag { "wf-humvar-${task.process}" }
//         cpus 4
//         memory '32 GB'
//         time '12h'
//         container params.epi2me_container ?: 'library/epi2me-wf:latest'
//         publishDir { params.publish_dir ?: "results/${project_name_ch}/${sample_name_ch}/wf-humvar" }, mode: 'copy'

//         input:
//             path bam_dir from bam_dir_ch
//             path ref from ref_ch
//             path targets from targets_ch
//             path tr from tr_ch
//             val sample_name from sample_name_ch
//             val project_name from project_name_ch

//         output:
//             path "results/${project_name}/${sample_name}/wf-humvar/*" into wf_outputs

//         script:
//             // --override_basecaller_cfg not included
//             """
//             mkdir -p results/${project_name}/${sample_name}/wf-humvar
//             # Run nested Nextflow pipeline using Singularity to provide Nextflow runtime
//             singularity exec --cleanenv ${container} \
//                 bash -lc "nextflow run epi2me-labs/wf-human-variation -r sha2b856c1f358ddf1576217a336bc0e9864b6dc0ed --bam '${bam_dir}' --ref '${ref}' --bed '${targets}' --out_dir results/${project_name}/${sample_name}/wf-humvar --sample_name '${sample_name}' --sv --snp --mod --str --phased --output_gene_summary --output_xam_fmt bam --modkit_args --preset 'traditional' --threads ${params.wf_humvar_threads ?: 8} -profile ${params.wf_humvar_profile ?: 'standard'}"
//             """
//     }

//     wf_humvar_out_ch = wf_outputs
// }
