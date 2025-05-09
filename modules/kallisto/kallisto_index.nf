process KALLISTO_INDEX {
    tag "$transcript_fasta"

    cpus 1
    memory 10.GB
    time '1:00:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/how-are-we-stranded-here:v1.0.1-e6ce74d'

    publishDir "${params.pubdir}/kallisto_index", pattern:"kallisto_index", mode:'copy'

    input:
        path(transcript_fasta)

    output:
        path("kallisto_index"), emit: kallisto_index

    script:
        """
        kallisto index -i kallisto_index --make-unique ${transcript_fasta}
        """
}

// Module from: https://github.com/KU-GDSC/workflows
