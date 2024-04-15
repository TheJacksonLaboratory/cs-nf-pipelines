process BCFTOOLS_VCF_TO_BCF {
    tag "$sampleID"

    cpus = 8
    memory = 6.GB
    time = '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), path(vcf), path(tbi)

    output:
    tuple val(sampleID), path("*bcf"), path("*csi"), emit: bcf

    script:

    """
    bcftools view ${vcf} -Ob -o ${vcf.baseName}.bcf
    bcftools index ${vcf.baseName}.bcf
    """
}
