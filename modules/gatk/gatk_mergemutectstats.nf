process GATK_MERGEMUTECTSTATS {
    tag "$sampleID"

    cpus = 1
    memory = 15.GB
    time '05:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    input:
    tuple val(sampleID), path(list)

    output:
    tuple val(sampleID), file("*.stats"), emit: stats

    script:
    //Estimate somatic variants using Mutect2
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    stats = list.collect { "-stats $it" }.join(' ')

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" MergeMutectStats \
    ${stats} \
    -O ${sampleID}_merged.stats
    """
}
