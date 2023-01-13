process RNA_SUMMARY_STATS {
    tag "$sampleID"

    cpus = 1
    time = '00:15:00'

    container 'quay.io/jaxcompsci/perl:0.1.0'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'summary_stats' }", pattern: "*stats.txt", mode:'copy'

    input:
    tuple val(sampleID), file(rsem_stats), file(quality_stats), file(picard_metrics)

    output:
    tuple val(sampleID), file("*.txt")

    script:
    log.info "----- Summary Metrics running on ${sampleID} -----"

    if (params.read_type == "PE")

      """
      perl ${projectDir}/bin/rnaseq/summary_QC_metrics_without_xenome.pl \
      ${quality_stats} \
      ${rsem_stats} \
      ${picard_metrics} > ${sampleID}_summary_stats.txt
      """

    else if (params.read_type == "SE")

      """
      perl ${projectDir}/bin/rnaseq/summary_QC_metrics_without_xenome_SE.pl \
      ${quality_stats} \
      ${rsem_stats} \
      ${picard_metrics}  > ${sampleID}_summary_stats.txt
      """
}
