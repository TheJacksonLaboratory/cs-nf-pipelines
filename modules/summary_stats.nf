process summ_stats {
    tag "sampleID"

    cpus = 1
    time = '00:15:00'
    clusterOptions = '-q batch'
    container

    publishDir "${sample_tmpdir}_tmp", pattern: "*stats.txt", mode: 'copy'

    input:
    tuple sampleID, file(fq_stat)
    tuple sampleID, file(aln_stat)
    tuple sampleID, file(mets_stat)

    output:
    tuple sampleID, file("*.txt")

    when:
    params.gen_org == "human"

    script:
    log.info "-----Human Summary Metrics running on ${sampleID}-----"

    if (params.reads == "PE")

      """
      perl ${params.summary_mets_PE} \
      ${fq_stat} \
      ${aln_stat} \
      ${mets_stat} > ${sampleID}_summary_stats.txt
      """

    else if (params.reads == "SE")

      """
      perl ${params.summary_mets_SE} \
      ${fq_stat} \
      ${aln_stat} \
      ${mets_stat}  > ${sampleID}_summary_stats.txt

      """

    }
