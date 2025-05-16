process REMOVE_DUPLICATE_READS {
    tag "$sampleID"

    cpus 2
    memory 4.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    input:
    tuple val(sampleID), file(marked_bam_file), file(marked_bai_file)

    output:
    tuple val(sampleID), file("*.sorted.rmDup.bam"), emit: rmDup_bam
    tuple val(sampleID), file("*.sorted.rmDup.bam.bai"), emit: rmDup_bai

    script:
    // Exclude reads flagged as pcr or optical duplicates (0x400), marked with bit flag 1024 in the BAM.
    """
    samtools view -h -b -F 1024 ${marked_bam_file} > ${sampleID}.sorted.rmDup.bam

    samtools index ${sampleID}.sorted.rmDup.bam
    """
}
