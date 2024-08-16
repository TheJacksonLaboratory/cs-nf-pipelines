process SAMTOOLS_FAIDX {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    publishDir "${params.pubdir}/g2gtools", pattern:"*.fai", mode:'copy', enabled: params.workflow == 'generate_pseudoreference' ? true : false
    publishDir "${params.pubdir}/genome_info", pattern:"*.fai", mode:'copy', enabled: params.workflow == 'chipseq' ? true : false

    input:
        tuple val(sampleID), path(fasta)

    output:
        tuple val(sampleID), path("*.fai"), emit: fai

    script:
    """
      samtools faidx ${fasta}
    """
}
