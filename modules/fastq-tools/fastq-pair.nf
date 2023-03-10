process FASTQ_PAIR {
    tag "$sampleID"

    cpus 1
    memory { 50.GB * task.attempt }
    time { 10.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 1

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
