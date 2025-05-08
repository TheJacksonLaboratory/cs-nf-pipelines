process PICARD_COLLECTALIGNMENTSUMMARYMETRICS{
    tag "$sampleID"

    cpus = 1
    memory = 20.GB
    time = '14:00:00'
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
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp"  CollectAlignmentSummaryMetrics \
    --INPUT ${bam} \
    --OUTPUT ${sampleID}_AlignmentMetrics.txt \
    --REFERENCE_SEQUENCE ${params.ref_fa} \
    --METRIC_ACCUMULATION_LEVEL ALL_READS \
    --VALIDATION_STRINGENCY LENIENT
    """
}
