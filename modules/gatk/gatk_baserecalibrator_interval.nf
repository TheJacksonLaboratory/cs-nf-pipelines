process GATK_BASERECALIBRATOR {
    tag "$sampleID"

    cpus = 1
    memory 30.GB
    time '08:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(bam), val(interval)

    output:
    tuple val(sampleID), path("*.table"), emit: table

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]
    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" BaseRecalibrator \
    -I ${bam} \
    -R ${params.ref_fa} \
    --known-sites ${params.dbSNP} \
    --known-sites ${params.gold_std_indels} \
    --known-sites ${params.phase1_1000G} \
    -L ${interval} \
    -O ${sampleID}_${interval}_recal_data.table
    """
}
