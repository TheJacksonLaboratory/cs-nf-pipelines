process FASTQ_SORT {
    tag "$sampleID"

    cpus 1
    memory { 50.GB * task.attempt }
    time { 2.h * task.attempt }
    errorStrategy 'retry'
    maxRetries 1

    container 'quay.io/biocontainers/fastq-tools:0.8.3--hbd632db_2'

    input:
    tuple val(sampleID), file(trimmed_hsa)

    output:
    tuple val(sampleID), file("*sorted_human*{1,2}.fastq"), emit: sorted_fastq

    script:
    command_two = params.read_type == 'PE' ? "fastq-sort --id ${trimmed_hsa[1]} > ${sampleID}_sorted_human_2.fastq" : ''

    """
    fastq-sort --id ${trimmed_hsa[0]} > ${sampleID}_sorted_human_1.fastq
    ${command_two}
    """
}
