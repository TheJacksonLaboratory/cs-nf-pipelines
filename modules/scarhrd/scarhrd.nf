process SCARHRD {
    tag "$sampleID"

    cpus 2
    memory 60.GB
    time '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/sequenza:v1'

    publishDir "${params.pubdir}/${sampleID}", pattern: "*HRD_score.txt", mode:'copy'


    input:
    tuple val(sampleID), val(meta), path(seqz)

    output:
    tuple val(sampleID), val(meta), path("*HRD_score.txt"), emit: hrd_score

    script:
    """
    Rscript ${projectDir}/bin/wes/scarHRD.R ${seqz} ${sampleID}
    """
}
