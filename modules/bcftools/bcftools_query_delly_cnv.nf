process BCFTOOLS_QUERY_DELLY_CNV {
    tag "$sampleID"
    
    cpus = 8
    memory = 6.GB
    time = '06:00:00'
    errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

    container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

    publishDir "${params.pubdir}/${sampleID + '/callers'}", pattern: "*.bed", mode: 'copy'

    input:
    tuple val(sampleID), path(bcf), path(csi), val(meta), val(normal_name), val(tumor_name), val(caller)

    output:
    tuple val(sampleID), path("*segmentation.bed"), val(meta), val(normal_name), val(tumor_name), val('delly_cnv'), emit: segmentation_file

    script:

    """
    bcftools query -s ${tumor_name} -H -f "%CHROM\\t%POS\\t%INFO/END\\t%ID\\t[%RDCN]\\t[%RDSD]\\t[%CN]\n" ${bcf} > ${sampleID}_delly_somatic_cnv_segmentation.bed
    """
}
