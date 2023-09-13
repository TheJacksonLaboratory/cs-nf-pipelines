process SURVIVOR_TO_BED {
    tag "$sampleID"

    cpus 1
    memory 100.GB
    time "00:30:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'rocker/tidyverse:4.2.1'

    input:
        tuple val(sampleID), file(annot), file(summary)
    output:
        tuple val(sampleID), file("${sampleID}.ins.bed"), file("${sampleID}.del.bed"), file("${sampleID}.dup.bed"), file("${sampleID}.inv.bed"), file("${sampleID}.tra.bed"), emit: sv_beds
    script:
        """
        /usr/bin/env Rscript ${projectDir}/bin/germline_sv/surv_annot_process.R ${sampleID}
        """
}
