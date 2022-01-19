process SUMMARY_STATS {
    tag "sampleID"

    cpus = 1
    time = '00:15:00'
    clusterOptions = '-q batch'
    
    container 'R_perl.sif'

    publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'summary_stats' }", pattern: "*stats.txt", mode:'copy'

    input:
    tuple val(sampleID), file(rsem_stats)
    tuple val(sampleID), file(quality_stats)
    tuple val(sampleID), file(picard_metrics)

    output:
    tuple val(sampleID), file("*.txt")

    when:
    params.gen_org == "human"

    script:
    log.info "----- Human Summary Metrics running on ${sampleID} -----"

    if (params.read_type == "PE")

      """
      perl ${params.summary_mets_PE} \
      ${quality_stats} \
      ${rsem_stats} \
      ${picard_metrics} > ${sampleID}_summary_stats.txt
      """

    else if (params.read_type == "SE")

      """
      perl ${params.summary_mets_SE} \
      ${quality_stats} \
      ${rsem_stat} \
      ${picard_metrics}  > ${sampleID}_summary_stats.txt

      """

    }
