process GRIDSS_CHROM_FILTER {
    tag "$sampleID"

    cpus = 1
    memory = 1.GB
    time = '01:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/jaxcompsci/internal_tools:v1.0'

    stageInMode = 'copy'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*_gridss_sv_unfiltered_chroms.vcf", mode:'copy', enabled: params.keep_intermediate

    input:
    tuple val(sampleID), path(vcf), val(meta), val(normal_name), val(tumor_name)
    val(chroms)

    output:
    tuple val(sampleID), path('*_gridss_sv_unfiltered_chroms.vcf'), val(meta), val(normal_name), val(tumor_name), emit: gridss_chrom_vcf
    
    script:
    chrom_list = chroms.collect { "$it" }.join(' ')

    """
    python ${projectDir}/bin/pta/filter_vcf.py \
    --vcf-file ${vcf} \
    --output ${sampleID}_gridss_sv_unfiltered_chroms.vcf \
    --chroms ${chrom_list}
    """
}
