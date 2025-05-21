process PBMM2_INDEX {

    tag "${fasta.baseName}"

    cpus 8
    memory 40.GB
    time "2:00:00"
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/pbmm2:1.9.0--h9ee0642_0'

    input:
        path(fasta)
    output:
        path "${fasta.baseName}.mmi", emit: pbmm2_index
    script:
        """
        pbmm2 index ${fasta} ${fasta.baseName}.mmi
        """
}
