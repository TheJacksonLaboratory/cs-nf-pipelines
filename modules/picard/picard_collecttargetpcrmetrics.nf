process PICARD_COLLECTTARGETPCRMETRICS {
    tag "$sampleID"

    cpus = 1
    memory = 5.GB
    time = '08:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.txt", mode:'copy'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.txt"), emit: txt

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" CollectTargetedPcrMetrics \
    --INPUT ${bam} \
    --OUTPUT ${sampleID}_CollectTargetedPcrMetrics.txt \
    --REFERENCE_SEQUENCE ${params.ref_fa} \
    --AMPLICON_INTERVALS ${params.amplicon_primer_intervals} \
    --TARGET_INTERVALS ${params.amplicon_target_intervals} \
    --COVERAGE_CAP 1500 \
    --NEAR_DISTANCE 50 \
    --VALIDATION_STRINGENCY LENIENT
    """
}
