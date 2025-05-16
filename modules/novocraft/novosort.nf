process NOVOSORT_MARKDUPLICATES {
    tag "$sampleID"

    cpus = 1
    memory = 8.GB
    time = '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/novosort:lastest'
    publishDir "${params.pubdir}/${sampleID}", pattern:"*_fixed_mate_dup_marked.bam", mode:'copy'

    input:
    tuple val(sampleID), file(fixed_mate_bam)

    output:
    tuple val(sampleID), file("*_fixed_mate_dup_marked.bam"), emit: fixed_mate_dup_marked_bam

    script:

    """
    novosort -markduplicates -t . -m 8G \
    ${fixed_mate_bam} > ${sampleID}_fixed_mate_dup_marked.bam
    """
}
