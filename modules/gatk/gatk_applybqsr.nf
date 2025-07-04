process GATK_APPLYBQSR {
    tag "$sampleID"

    cpus = 1
    memory = 60.GB
    time = '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/bam'}", pattern: "*.ba*", mode:'copy'

    input:
    tuple val(sampleID), file(bam), file(table)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam
    tuple val(sampleID), file("*.bai"), emit: bai

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" ApplyBQSR \
    -R ${params.ref_fa} \
    -I ${bam} \
    --bqsr-recal-file ${table} \
    -O ${sampleID}_realigned_BQSR.bam
    """
}
