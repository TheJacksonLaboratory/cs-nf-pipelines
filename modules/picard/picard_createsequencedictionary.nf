process PICARD_CREATESEQUENCEDICTIONARY {
    tag "${fasta}"

    cpus 1
    memory 24.GB
    time '06:00:00'

    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container "quay.io/biocontainers/picard:2.26.10--hdfd78af_0"

    publishDir "${params.pubdir}", pattern: "${fasta.baseName}.dict", mode:'copy'

    input:
        path(fasta)

    output:
        path("${fasta.baseName}.dict"), emit: dict

    script:
        String my_mem = (task.memory-1.GB).toString()
        my_mem =  my_mem[0..-4]

        """
        picard -Xmx${my_mem}G CreateSequenceDictionary \
        REFERENCE=${fasta} \
        OUTPUT=${fasta.baseName}.dict
        """
    }

// Module from: https://github.com/KU-GDSC/workflows
