process FILTER_REMOVE_MULTI_SIEVE {
    tag "$sampleID"

    cpus 8
    memory 20.GB
    time '04:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

    input:
    tuple val(sampleID), file(shift_bams)

    output:
    tuple val(sampleID), file("*.shift.tmp.ba*")

    script:
    """
    alignmentSieve \
    --numberOfProcessors $task.cpus \
    --ATACshift \
    --bam ${shift_bams[0]} \
    -o ${sampleID}.shift.tmp.bam 
    """
}
