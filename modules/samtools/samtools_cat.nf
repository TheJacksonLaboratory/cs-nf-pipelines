process SAMTOOLS_CAT {
    tag "$sampleID"

    cpus 1
    memory 8.GB
    time '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

    input:
    tuple val(sampleID), path(bams)

    output:
    tuple val(sampleID), path("*.bam"), emit: bam

    script:

        """
        echo ${bams} > FOFN.txt
        samtools cat -b FOFN.txt -o ${sampleID}_concatBAM.bam
        """
}
