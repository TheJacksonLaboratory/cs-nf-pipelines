process SAMTOOLS_MPILEUP {
    tag "$sampleID"

    cpus 2
    memory 60.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    input:
    tuple val(sampleID), val(meta), path(bam), path(bai)

    output:
    tuple val(sampleID), val(meta), file("*.pileup.gz"), emit: pileup

    script:
    """
    samtools mpileup -f ${params.ref_fa} -Q 20 ${bam} | gzip -c - > ${bam.baseName}.pileup.gz
    """
}
