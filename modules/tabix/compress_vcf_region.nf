process COMPRESS_INDEX_VCF_REGION {
    tag "$sampleID"

    cpus = 1
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/tabix:1.11--hdfd78af_0'

    input:
    tuple val(sampleID), file(vcf), val(meta), val(chrom)

    output:
    tuple val(sampleID), file("*.vcf.gz"), file("*.vcf.gz.tbi"), val(meta), val('empty_name'), val('empty_name'), val(chrom), emit: compressed_vcf_tbi

    """
    bgzip \
    -c \
    ${vcf} \
    > ${vcf}.gz

    tabix ${vcf}.gz
    """
}
