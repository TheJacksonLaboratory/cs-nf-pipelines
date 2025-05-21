process BCFTOOLS_SPLITMULTIALLELIC {
    tag "$sampleID"

    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), file(vcf), file(tbi), val(meta), val(normal_name), val(tumor_name), val(tool)

    output:
    tuple val(sampleID), file("*.vcf"), val(meta), val(normal_name), val(tumor_name), val(tool), emit: vcf

    script:
    output_name = vcf.getBaseName().replace('.vcf', '')
    """
    bcftools \
    norm \
    -m \
    -any \
    --threads ${task.cpus} \
    --no-version \
    -f ${params.ref_fa} \
    -o ${output_name}_multiAllelicSplit.vcf \
    ${vcf}
    """
}
