process XENGSORT_INDEX {

    // resource utilization
    cpus 32
    memory 60.GB
    time 1.h
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    // load xengsort container
    container 'quay.io/biocontainers/xengsort:2.0.5--pyhdfd78af_0'

    // output directory
    publishDir "${params.pubdir}/xengsort/", mode: 'copy'

    // inputs
    input:
    path(xengsort_host_fasta)
    path(graft_fasta)

    output:
    // index output
    path("xengsort_index/"), emit: xengsort_index

    script:
    """
    xengsort index \
    --index ${params.xengsort_idx_name} \
    -H ${xengsort_host_fasta} \
    -G ${graft_fasta} \
    -k 25 \
    -n 4_500_000_000 \
    -W ${task.cpus}

    mkdir xengsort_index
    mv *.info *.hash xengsort_index/
    """
}
