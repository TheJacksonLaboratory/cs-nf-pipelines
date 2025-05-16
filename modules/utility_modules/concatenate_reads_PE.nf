process CONCATENATE_READS_PE {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'ubuntu:20.04'

    publishDir "${params.pubdir}/${sampleID + '/concatenated_reads'}", pattern: "*", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), file(R1), file(R2)

    output:
    tuple val(sampleID), file("*"), emit: concat_fastq

    script:

    """
    cat $R1 > ${sampleID}_R1${params.extension}
    cat $R2 > ${sampleID}_R2${params.extension}
    """
}
