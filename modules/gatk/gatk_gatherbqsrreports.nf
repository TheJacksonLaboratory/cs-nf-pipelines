process GATK_GATHERBQSRREPORTS {
    tag "$sampleID"

    cpus = 1
    memory = 40.GB
    time = '03:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'broadinstitute/gatk:4.2.4.1'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern: "*.table", mode:'copy'

    input:
    tuple val(sampleID), path(reports)

    output:
    tuple val(sampleID), path("*.table"), emit: table

    script:
    String my_mem = (task.memory-1.GB).toString()
    my_mem =  my_mem[0..-4]

    inputs = reports.collect { "-I $it" }.join(' ')

    """
    mkdir -p tmp
    gatk --java-options "-Xmx${my_mem}G -Djava.io.tmpdir=`pwd`/tmp" GatherBQSRReports \
    ${inputs} \
    -O ${sampleID}_recal_data.table
    """
}
