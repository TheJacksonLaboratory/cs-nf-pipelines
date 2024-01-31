process BEDOPS_SORT {
    tag "Sorting Bed"

    cpus 1
    memory 25.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bedops:2.4.41--h4ac6f70_1'

    input:
    path(bed)

    output:
    path("*sorted.bed"), emit: sorted_bed
 
    script:
    """
    sort-bed ${bed} > ${bed.baseName}.sorted.bed
    """
}
