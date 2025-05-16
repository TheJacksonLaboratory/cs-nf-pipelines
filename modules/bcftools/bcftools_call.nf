process BCFTOOLS_CALL {
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

    bcftools call \
    --output-type v \
    --output ${sampleID}.mpileup.called.vcf \
    --multiallelic-caller \
    ${vcf}.gz
    """
}
