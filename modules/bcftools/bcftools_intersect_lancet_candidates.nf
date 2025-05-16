process BCFTOOLS_INTERSECTVCFS {
    tag "$sampleID"

    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), file(candidate_vcf), file(candidate_tbi), file(lancet_confirm_vcf), file(lancet_confirm_tbi), val(meta), val(normal_name), val(tumor_name), val(chrom)

    output:
    tuple val(sampleID), file("*.vcf.gz"), file("*.tbi"), val(meta), val(normal_name), val(tumor_name), emit: vcf

    script:
    """
    bcftools \
    isec \
    -w 1 \
    -c none \
    -n =2 \
    --threads ${task.cpus} \
    ${lancet_confirm_vcf} \
    ${candidate_vcf} \
    > ${sampleID}_confirmed_lancet_merged_${chrom}.vcf

    bgzip ${sampleID}_confirmed_lancet_merged_${chrom}.vcf
    tabix ${sampleID}_confirmed_lancet_merged_${chrom}.vcf.gz
    """
}
