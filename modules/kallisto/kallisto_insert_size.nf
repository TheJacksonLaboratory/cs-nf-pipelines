process KALLISTO_INSERT_SIZE {
    tag "$sampleID"

    cpus 1
    memory 1.GB
    time '00:05:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    cache 'lenient'

    container 'quay.io/biocontainers/pizzly:0.37.3--h470a237_3'

    input:
        tuple val(sampleID), val(kallisto_abundance)

    output:
        tuple val(sampleID), path('insert_size.txt'), emit: kallisto_insert_size

    script:
    """
    python ${projectDir}/bin/rna_fusion/compute_insert_size.py ${kallisto_abundance} > insert_size.txt
    """
}
