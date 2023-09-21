    process SMOOVE_CALL {
    tag "$sampleID"

    cpus = 12
    memory 80.GB
    time '18:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.mem} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'brentp/smoove:v0.2.7'

    input:
    tuple val(sampleID), file(bam_bwa_lumpy_sort), file(bam_bwa_lumpy_sort_bai)

    output:
    tuple val(sampleID), file("${sampleID}*.vcf.gz"), emit: lumpy_vcf

    script:
    """
    export TMPDIR=./

    smoove call \
        -x \
        --name ${sampleID} \
        --exclude ${params.exclude_regions} \
        --fasta ${params.fasta} \
        --processes ${task.cpus} \
        --support 3 \
        --genotype ${bam_bwa_lumpy_sort}

    """
}
