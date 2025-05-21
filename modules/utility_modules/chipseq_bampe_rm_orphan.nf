process BAMPE_RM_ORPHAN {
    tag "$sampleID"

    cpus 1
    memory 15.GB
    time '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/mulled-v2-57736af1eb98c01010848572c9fec9fff6ffaafd:402e865b8f6af2f3e58c6fc8d57127ff0144b2c7-0'

    input:
    tuple val(sampleID), file(bam)

    output:
    tuple val(sampleID), file("*.bam"), emit: bam

    script:  // This script was bundled withing the nf-core/chipseq/bin/ directory
    prefix = "${sampleID}.mLb.clN"
    """
    python ${projectDir}/bin/chipseq/bampe_rm_orphan.py ${bam[0]} ${prefix}.bam --only_fr_pairs
    """
}
