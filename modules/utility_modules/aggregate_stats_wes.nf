process AGGREGATE_STATS {
  tag "$sampleID"

  cpus = 1
  time = '00:30:00'
  clusterOptions = '-q batch'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'aggregate_stats' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(filter_stats)
  tuple val(sampleID), file(algn_met)
  tuple val(sampleID), file(picard_met)

  output:
  tuple val(sampleID), file("*summary_stats.txt"), emit: txt

  script:
  log.info "----- Generating Summary Stats for: ${sampleID} -----"

  """
  python ${projectDir}/bin/wes/aggregate_stats_wes.py ${sampleID}_summary_stats.txt ${filter_stats} ${picard_met} ${algn_met}
  """
}
