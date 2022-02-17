process AGGREGATE_STATS {
  tag "sampleID"

  cpus = 1
  time = '00:30:00'
  clusterOptions = '-q batch'

  container 'python_2.7.3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'aggregate_stats' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(filter_stats)
  tuple val(sampleID), file(algn_met)
  tuple val(sampleID), file(picard_met)

  output:
  tuple val(sampleID), file("*summary_stats.txt"), emit: txt

  script:
  log.info "----- Generating Summary Stats for: ${sampleID} -----"

  """
  python ${params.stats_agg} ${sampleID}_summary_stats.txt ${filter_stats} ${picard_met} ${algn_met}
  """
}
