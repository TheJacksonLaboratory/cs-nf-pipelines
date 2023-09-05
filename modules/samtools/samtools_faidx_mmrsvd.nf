process SAMTOOLS_FAIDX {
    tag "${fasta}"

    cpus 1
    memory 4.GB
    time '2:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/samtools:1.10--h9402c20_2'

    input:
        file(fasta)

    output:
        tuple file("${fasta}"), file("${fasta}.fai"), emit: fasta_fai

    script:
        """
        samtools faidx ${fasta}
        """
}