process SEQUENZA_NA_WINDOWS {
    tag "$sampleID"

    cpus = 1
    memory = 5.GB
    time = '00:45:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/sequenza:v1'

    input:
    tuple val(sampleID), val(meta), path(extract_rdata)

    output:
    tuple val(sampleID), path("*win.txt"), emit: na_windows

    script:
    """
    Rscript ${projectDir}/bin/wes/sequenza_seg_na_window.R ${extract_rdata} ${sampleID}_win.txt
    """
}
