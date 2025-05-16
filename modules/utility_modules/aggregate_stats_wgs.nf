process AGGREGATE_STATS {
    tag "$sampleID"

    cpus = 1
    time = '00:30:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

    publishDir "${params.pubdir}/${sampleID + '/stats'}", pattern:"*.txt", mode:'copy'

    input:
    tuple val(sampleID), file(filter_stats), file(picard_met), file(algn_met), file(cov_met)

    output:
    tuple val(sampleID), file("*summary_stats.txt"), emit: txt

    script:

    """
    python ${projectDir}/bin/wgs/aggregate_stats_wgs.py ${sampleID}_summary_stats.txt ${filter_stats} ${picard_met} ${algn_met} ${cov_met}
    """
}
