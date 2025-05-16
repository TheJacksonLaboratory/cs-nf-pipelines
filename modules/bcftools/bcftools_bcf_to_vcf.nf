process BCFTOOLS_BCF_TO_VCF {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.vcf.gz", mode: 'copy'

    input:
    tuple val(sampleID), path(bcf), path(csi), val(meta), val(normal_name), val(tumor_name), val(caller)

    output:
    tuple val(sampleID), path("*vcf.gz"), path("*tbi"), val(meta), val(normal_name), val(tumor_name), val(caller), emit: vcf_tbi

    script:

    """
    bcftools view ${bcf} -Oz -o ${bcf.baseName}.vcf.gz
    bcftools index -t ${bcf.baseName}.vcf.gz
    """
}
