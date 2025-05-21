process CUTADAPT {
    tag "$sampleID"

    cpus 8
    memory 10.GB
    time '20:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.log", mode: 'copy'

    container 'quay.io/biocontainers/cutadapt:2.3--py37h14c3975_0'

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*paired_trimmed.fq"), emit: paired_trimmed_fastq
    tuple val(sampleID), path("*.log"), emit: cutadapt_log

    script:

    paired_end = params.read_type == 'PE' ?  "-p ${sampleID}_R2_paired_trimmed.fq" : ''

    """
    cutadapt \
    -a ${params.cutadaptAdapterR1} \
    -A ${params.cutadaptAdapterR2} \
    --minimum-length ${params.cutadaptMinLength} \
    --quality-cutoff ${params.cutadaptQualCutoff} \
    -j $task.cpus \
    -o ${sampleID}_R1_paired_trimmed.fq \
    ${paired_end} \
    ${fq_reads} \
    > ${sampleID}_cutadapt.log \
    2>&1
    """
}
