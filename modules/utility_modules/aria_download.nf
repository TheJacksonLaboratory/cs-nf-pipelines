process ARIA_DOWNLOAD {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '10:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/aria2:1.36.0'

    input:
    tuple val(sampleID), val(meta), val(read_num), val(link)

    output:
    tuple val(sampleID), val(meta), val(read_num), path("*"), emit: file

    script:

    """
    aria2c --connect-timeout=180 --retry-wait=60 --timeout=180 ${link}
    """
}
