process BWA_INDEX {
    tag "${fasta}"

    cpus 1
    memory 30.GB
    time '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}
    
    container 'quay.io/biocontainers/bwa:0.7.17--hed695b0_6'

    input:
        file(fasta)
    output:
        path("${fasta}.*"), emit: bwa_index

    script:
        """
        bwa index ${fasta}
        """
}