process BCFTOOLS_DUPHOLD_FILTER {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '02:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    publishDir "${params.pubdir}/${sampleID + '/duphold'}", pattern:"*.vcf", mode:'copy'

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), path(vcf)
    val(tool)

    output:
    tuple val(sampleID), path("*vcf"), emit: vcf

    script:

    """
    bcftools view -e '(SVTYPE = "DEL" & FMT/DHFFC[0] > 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] < 1.3)' $vcf > ${sampleID}_${tool}_dupholdFiltered.vcf
    """
}
