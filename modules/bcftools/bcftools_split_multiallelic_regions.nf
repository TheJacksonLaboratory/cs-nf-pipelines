process BCFTOOLS_SPLITMULTIALLELIC_REGIONS {
    tag "$sampleID"

    cpus = 4
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    input:
    tuple val(sampleID), file(vcf), file(index)
    val(chrom_list)

    output:
    tuple val(sampleID), file("*.vcf.gz"), file("*.vcf.gz.tbi"), emit: vcf_idx

    script:
    
    listOfChroms = chrom_list.collect { "$it" }.join(',')

    """
    bcftools \
        norm \
        -m \
        -any \
        --threads ${task.cpus} \
        --regions ${listOfChroms} \
        --no-version \
        -f ${params.ref_fa} \
        -o ${sampleID}_split.vcf \
        ${vcf}

    bgzip \
        -c \
        ${sampleID}_split.vcf > ${sampleID}_split.vcf.gz

    tabix ${sampleID}_split.vcf.gz

    """
}
