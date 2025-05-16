process CONCATENATE_READS_SAMPLESHEET {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/concatenated_reads'}", pattern: "*fastq.gz", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), val(num_lanes), val(meta), val(read_num), path(reads)

    output:
    tuple val(sampleID), val(num_lanes), val(meta), val(read_num), path("*fastq.gz"), emit: concat_fastq

    when:
    num_lanes > 1

    script:

    """
    cat $reads > ${sampleID}_${read_num}.fastq.gz
    """
}
