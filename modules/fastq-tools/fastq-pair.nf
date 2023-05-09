process FASTQ_PAIR {
    tag "$sampleID"

    cpus 1
    memory 50.GB
    time { reads[0].size() < 35.GB ? 10.h : 18.h }
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/fastq-pair:1.0--h87f3376_3'

    input:
    tuple val(sampleID), file(reads)

    output:
    tuple val(sampleID), file("*.paired.fq"), emit: paired_fastq
    tuple val(sampleID), file("*.single.fq"), emit: single_fastq

    script:

    """
    fastq_pair ${reads[0]} ${reads[1]}
    """
}
