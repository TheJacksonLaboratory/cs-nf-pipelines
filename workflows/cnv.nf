#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Import modules
include { IAAP_CLI } from "${projectDir}/modules/illumina/iaap_cli"
include { help } from "${projectDir}/bin/help/cnv.nf"
include { param_log } from "${projectDir}/bin/log/cnv.nf"

// Help if needed
if (params.help) {
    help()
    exit 0
}
// Log parameter info
param_log()
// Parameter validation
if (!params.idat_folder || !params.bpm_file || !params.egt_file) {
    exit 1, "All parameters (idat_folder, bpm_file, egt_file) are required."
}
errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
// Main workflow
workflow CNV_ARRAY {
    IAAP_CLI(
        idat_folder: params.idat_folder
    )
}