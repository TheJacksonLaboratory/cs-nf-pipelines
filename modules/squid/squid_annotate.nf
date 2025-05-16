process SQUID_ANNOTATE {
    tag "$sampleID"

    cpus 1
    memory 10.GB
    time 5.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a'

    publishDir "${params.pubdir}/${sampleID + '/fusions'}", pattern: "*.{tsv,txt}", mode:'copy'

    input:
        tuple val(sampleID), path(txt)
        path(gtf)

    output:
        tuple val(sampleID), path("*annotated.txt"), emit: squid_fusions_annotated

    script:
    """
    AnnotateSQUIDOutput.py ${gtf} ${txt} ${sampleID}_squid_fusions_annotated.txt
    """
}
