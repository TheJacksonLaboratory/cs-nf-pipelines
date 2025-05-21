    process SMOOVE_CALL {
    tag "$sampleID"

    cpus = 12
    memory 80.GB
    time '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'brentp/smoove:v0.2.7'
    
    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.vcf.gz", mode: 'copy'

    input:
    tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name)

    output:
    tuple val(sampleID), path("*vcf.gz"), path("*vcf.gz.csi"), val(meta), val(normal_name), val(tumor_name), val('lumpy'), emit: vcf_tbi

    script:
    """
    export TMPDIR=./

    smoove call \
        -x \
        --name ${sampleID} \
        --exclude ${params.exclude_list} \
        --fasta ${params.ref_fa} \
        --processes ${task.cpus} \
        --support 5 \
        --duphold \
        --genotype ${tumor_bam} ${normal_bam}
    """
}

// NOTE: support: 5 was set following testing. Support 3 produced a large number of spurious calls in a t/n context. 
