process BCFTOOLS_ANNOTATE {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), path(vcf)

    output:
    tuple val(sampleID), path("*.vcf"), emit: vcf

    script:
    """
    bgzip -c -f ${vcf} > ${vcf}.gz

    bcftools index -f -t ${vcf}.gz

    bcftools annotate \
    --annotations ${params.snp_annotations} \
    --columns CHROM,FROM,TO,ID \
    ${vcf}.gz | bcftools view \
    -i "FILTER='PASS' && ID=@${params.snpID_list}" \
    --output-file ${sampleID}.mpileup.called.filtered.annotated.vcf
    """
}
