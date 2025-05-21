process FILTER_REMOVE_MULTI_SHIFT {
    tag "$sampleID"

    cpus 4
    memory 10.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.sorted.rmDup.rmChrM.rmMulti.filtered.ba*", mode: 'copy' 

    input:
    tuple val(sampleID), file(mtdna_bam_file), file(mtdna_bai_file)

    output:
    tuple val(sampleID), file("*.shift.tmp0.ba*")
    tuple val(sampleID), file("*.sorted.rmDup.rmChrM.rmMulti.filtered.ba*"), emit: srf_bam

    script:
    // Filter reads unmapped, mate unmapped, not primary alignment, reads failing platform, pcr duplicates (-F 1804) and reatin properly paired reads (-f 2) in bam file
    """
    # filter low quality reads
    samtools view -@ $task.cpus -h -q 30 ${mtdna_bam_file} \
    > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.bam

    # filter reads unmapped, mate unmapped, not primary alignment, reads failing platform, pcr duplicates (-F 1804)
    # retain properly paired reads (-f 2)
    samtools view -@ $task.cpus -h -b -F 1804 -f 2 \
    ${sampleID}.sorted.rmDup.rmChrM.rmMulti.bam \
    > ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam

    samtools index \
    ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam

    samtools sort \
    -@ $task.cpus -O bam \
    -o ${sampleID}.shift.tmp0.bam \
    ${sampleID}.sorted.rmDup.rmChrM.rmMulti.filtered.bam

    samtools index \
    ${sampleID}.shift.tmp0.bam
    """
}
