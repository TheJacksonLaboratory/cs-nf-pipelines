process NON_CHAIN_REINDEX {
    tag "$sampleID"

    cpus  1
    memory 8.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.filtered.shifted.*", mode: 'copy'
    
    input:
    tuple val(sampleID), file(bam_shifted)

    output:
    tuple val(sampleID), file("*.filtered.shifted.*")

    when: params.chain == null

    script:
    // This module is for Reference Strain Samples.
    // To filter Mitochondrial, Unplaced/Unlocalized Reads from bam file.
    """
    samtools index ${bam_shifted[0]}

    # filter Mitochondrial, Unplaced/Unlocalized Reads
    samtools view ${bam_shifted[0]} \
    -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    | grep -ve 'SN:MT*\\|SN:GL*\\|SN:JH:*' > tmp.sam

    # Re-sort BAM following MT removal. This is required for PICARD MERGE
    samtools sort -@ $task.cpus tmp.sam \
    > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10.bam

    # Index BAM
    samtools index \
    ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10.bam
    """
}
