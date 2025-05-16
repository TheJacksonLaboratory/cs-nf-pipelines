process FILTER_GTF {
    tag "Filtering GTF"

    cpus 1
    memory 5.GB
    time 2.hour
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/ripgrep:v1'

    publishDir "${params.pubdir}", pattern: '*.gtf', mode:'copy'

    output:
    path('*.gtf'), emit: filtered_gtf

    script:
    include_statement = params.gtf_biotype_include.split(',').collect { "$it" }.join('|')
    """
    sh ${projectDir}/bin/g2gtools/filter_gtf.sh ${params.primary_reference_gtf} "${include_statement}"
    """

    stub:
    """
    touch "gtf.included.gtf"
    """
}
