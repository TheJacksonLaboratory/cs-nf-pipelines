process ASCAT {
    tag "$sampleID"

    cpus 1
    memory 24.GB
    time '01:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/ascat:v3.1.3'

    publishDir "${params.pubdir}/${sampleID}", mode: 'copy'

    input:
        tuple val(sampleID), val(meta), path(BAF), path(LRR)

    output:
        tuple val(sampleID), val(meta), path("*.txt"), emit: all_txt
        tuple val(sampleID), val(meta), path("*.png"), emit: all_png
        tuple val(sampleID), val(meta), path("*.Rdata"), emit: ascat_rdata
        tuple val(sampleID), val(meta), path("*segments_raw.txt"), path("*.ploidy.txt"), emit: seg_ploidy

    script:
        """
        Rscript ${projectDir}/bin/cnv_array/ASCAT_run.R \
            ${sampleID} \
            ${BAF} \
            ${LRR} \
            ${meta.gender} \
            ${params.snp_platform} \
            ${params.gc_file} \
            ${params.rt_file}
        """
}
