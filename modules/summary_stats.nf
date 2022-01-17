process SUMMARY_STATS {
    tag "sampleID"

    cpus = 1
    time = '00:15:00'
    clusterOptions = '-q batch'
    container

    publishDir "${pubdir}/summary_stats", pattern: "*stats.txt", mode: 'copy'

    input:
    tuple sampleID, file(rsem_stats)
    tuple sampleID, file(quality_stats)
    tuple sampleID, file(picard_metrics)

    output:
    tuple sampleID, file("*.txt")

    when:
    params.gen_org == "human"

    script:
    log.info "-----Human Summary Metrics running on ${sampleID}-----"

    if (params.reads == "PE")

      """
      perl ${params.summary_mets_PE} \
      ${quality_stats} \
      ${aln_stat} \
      ${mets_stat} > ${sampleID}_summary_stats.txt
      """

    else if (params.reads == "SE")

      """
      perl ${params.summary_mets_SE} \
      ${quality_stats} \
      ${aln_stat} \
      ${mets_stat}  > ${sampleID}_summary_stats.txt

      """

    }
