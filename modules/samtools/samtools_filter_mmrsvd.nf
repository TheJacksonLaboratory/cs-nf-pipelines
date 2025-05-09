process SAMTOOLS_FILTER {
    tag "${sampleID}"

    cpus 1
    memory 50.GB
    time '10:00:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    publishDir "${params.pubdir}/${sampleID + '/alignments'}", mode:'copy'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("${sampleID}.q30.bam"), file("${sampleID}.q30.bam.bai"), emit: bam_and_index

    script:
        """
        samtools view -bq 30 ${bam} -o ${sampleID}.q30.bam
        samtools index ${sampleID}.q30.bam
        """
}
