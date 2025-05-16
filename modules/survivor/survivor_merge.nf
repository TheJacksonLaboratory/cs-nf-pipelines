process SURVIVOR_MERGE {
    tag "$sampleID"

    cpus 8
    memory 40.GB
    time "4:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/survivor:1.0.7--he513fc3_0'

    input:
        tuple val(sampleID), file(vcf_tuple)
    output:
        tuple val(sampleID), file("${sampleID}_mergedCall.*.vcf"), emit: vcf
    script:
        if (params.data_type == "pacbio")
            """
            ls ${vcf_tuple[0]} > vcf_list.txt
            ls ${vcf_tuple[1]} >> vcf_list.txt
            SURVIVOR merge vcf_list.txt ${params.surv_dist} ${params.surv_supp} ${params.surv_type} ${params.surv_strand} 0 ${params.surv_min} ${sampleID}_mergedCall.PS.vcf
            """
        else if (params.data_type == "illumina")
            """
            ls ${vcf_tuple[0]} > vcf_list.txt
            ls ${vcf_tuple[1]} >> vcf_list.txt
            ls ${vcf_tuple[2]} >> vcf_list.txt
            SURVIVOR merge vcf_list.txt ${params.surv_dist} ${params.surv_supp} ${params.surv_type} ${params.surv_strand} 0 ${params.surv_min} ${sampleID}_mergedCall.DLM.vcf
            """
        else if (params.data_type == "ont")
            """
            ls ${vcf_tuple[0]} > vcf_list.txt
            ls ${vcf_tuple[1]} >> vcf_list.txt
            SURVIVOR merge vcf_list.txt ${params.surv_dist} ${params.surv_supp} ${params.surv_type} ${params.surv_strand} 0 ${params.surv_min} ${sampleID}_mergedCall.NS.vcf
            """            
}
