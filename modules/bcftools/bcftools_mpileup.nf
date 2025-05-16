process BCFTOOLS_MPILEUP {
    tag "$sampleID"
    
    cpus = 8
    memory = 20.GB
    time = '12:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("*.vcf"), emit: vcf

    script:
    regions_file = params.genotype_targets ? "--regions-file ${params.genotype_targets}" : ""
    """
    bcftools mpileup \
    --output ${sampleID}.mpileup.vcf \
    --fasta-ref ${params.ref_fa} \
    --min-MQ 30 \
    --min-BQ 20 \
    ${regions_file} \
    --skip-indels \
    ${bam}
    """
}
