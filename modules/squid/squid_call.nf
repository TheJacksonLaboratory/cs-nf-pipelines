process SQUID {
    tag "$sampleID"

    cpus 1
    memory 50.GB
    time 5.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}


    container 'docker.io/nfcore/rnafusion:squid_1.5-star2.7.1a'

    input:
        tuple val(sampleID), path(bam), path(chimeric_bam)

    output:
        tuple val(sampleID), path("*sv.txt"), emit: squid_fusions

    script:
    """
    squid -b ${bam} -c ${chimeric_bam} -o ${sampleID}.squid.fusions
    """
}
