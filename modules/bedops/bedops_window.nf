process BEDOPS_WINDOW {
    tag "Window Bed"

    cpus 1
    memory 25.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bedops:2.4.41--h4ac6f70_1'

    input:
    path(sorted_bed)
    path(windows)

    output:
    path("*WindowCount.bed"), emit: window_bed
 
    script:
    """
    bedops -m ${sorted_bed} > Sorted_Whole_Exome_Capture.merged.bed; bedmap --echo --bases --delim '\t' ${windows} Sorted_Whole_Exome_Capture.merged.bed > ${windows.baseName}.${sorted_bed.baseName}.WindowCount.bed
    """
}
