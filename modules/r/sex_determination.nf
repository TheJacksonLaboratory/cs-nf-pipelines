process SEX_DETERMINATION {
    tag "$sampleID"

    cpus 1
    memory 5.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'rocker/r-ver:4.4.1'

    publishDir "${params.pubdir}/${sampleID}", mode:'copy'

    input:
        tuple val(sampleID), path(counts)
    output:
        tuple val(sampleID), file("*.csv"), emit: csv

    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/rnaseq/sex_determination.R ${counts} ./ ${sampleID}
        """
}
