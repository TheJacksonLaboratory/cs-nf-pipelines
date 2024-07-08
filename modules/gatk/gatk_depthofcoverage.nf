process GATK_DEPTHOFCOVERAGE {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), file(bam), file(bai)
    val(L)

    output:
    tuple val(sampleID), file("*_gatk_temp.txt"), emit: txt

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" DepthOfCoverage \
    -R ${params.ref_fa} \
    --output-format TABLE \
    -O ${sampleID}_gatk_temp.txt \
    -I ${bam} \
    -L  ${L} \
    --omit-per-sample-statistics \
    --omit-interval-statistics \
    --omit-locus-table \
    """
}
