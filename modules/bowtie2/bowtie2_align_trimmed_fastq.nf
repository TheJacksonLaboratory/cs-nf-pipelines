process ALIGN_TRIMMED_FASTQ {
    tag "$sampleID"

    cpus 16
    memory 30.GB
    time '48:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.log", mode: 'copy'
    container 'biocontainers/bowtie2:v2.4.1_cv1'

    input:
    tuple val(sampleID), path(fq_reads)

    output:
    tuple val(sampleID), path("*.sam"), emit: sam
    tuple val(sampleID), path("*_bowtie2.log"), emit: bowtie_log

    script:
    String options = params.bowtieVSensitive  == 'true' ? '--very-sensitive' : ''
    """
    bowtie2 \
    ${options} \
    -X ${params.bowtieMaxInsert} \
    -q \
    -p $task.cpus \
    -x ${params.bowtie2Index} \
    -1 ${fq_reads[0]} \
    -2 ${fq_reads[1]} \
    -S ${sampleID}.sam \
    2>${sampleID}_bowtie2.log \
    """
}
