process CHAIN_SORT_FIXMATE_BAM {
    tag "$sampleID"

    cpus  8
    memory 20.GB
    time '20:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.filtered.shifted.*", mode: 'copy'
    
    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), path("*.filtered.shifted.*")

    when: params.chain != null

    script:
    // This module is for Non-Reference Strain Samples. 
    // To sort bam by read name, fix the mate information, re-sort by coordinates and filter Mitochondrial Reads from bam file. 
    """
    # sort bam by read name
    samtools sort \
    -n \
    -@ $task.cpus -O bam \
    -o ${sampleID}.tmp3.mm10.bam ${bam[0]}

    # fix the mate information. This is done to fix 'TLEN' which is required for MACS2
    samtools fixmate \
    -O bam ${sampleID}.tmp3.mm10.bam ${sampleID}.tmp4.mm10.bam

    # re-sort bam by coordinates. This step is required for the MT filter to work properly
    samtools sort \
    -@ $task.cpus -O bam \
    -o ${sampleID}.tmp5.mm10.bam ${sampleID}.tmp4.mm10.bam

    samtools index ${sampleID}.tmp5.mm10.bam

    # filter Mitochondrial Reads
    samtools view ${sampleID}.tmp5.mm10.bam \
    -h 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y \
    | grep -ve 'SN:MT*' > tmp.sam

    # re-sort following MT removal. This is required for PICARD MERGE
    samtools sort -@ $task.cpus tmp.sam \
    > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10.bam

    # Index BAM
    samtools index \
    ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.shifted.mm10.bam
    """
}
