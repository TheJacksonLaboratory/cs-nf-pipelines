process AGGREGATE_STATS { // NEED TO FIX THIS UP TO TAKE N NUMBER OF FILES OR A LIST OF FILES
  tag "sampleID"

  cpus = 1
  time = '00:30:00'
  clusterOptions = '-q batch'

  container 'python_2.7.3.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(filter_stats), file(AlignMet) from fq_stats.join(covmet)

  output:
  tuple sampleID, file("*summary_stats.txt")
  tuple sampleID, file("*stats.txt") into dummy_stats_file

  script:
  log.info "-----Generating summary stats file for ${sampleID}-----"

  """
  python ${params.stats_agg} ${sampleID}_summary_stats.txt ${filter_stats} ${AlignMet}
  """
}
process AGGREGATE_STATS_HUMAN { // THIS IS DEPRECATED
  tag "sampleID"

  cpus = 1
  time = '00:30:00'
  clusterOptions = '-q batch'

  container 'python_2.7.3.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*.*", mode: 'copy'

  input:
  tuple sampleID, file(filter_stats), file(PicardMet), file(CoverageMet) from fq_stats.join(picard_metrics).join(covmet)

  output:
  tuple sampleID, file("*summary_stats.txt")
  tuple sampleID, file("*stats.txt") into dummy_stats_file

  script:
  log.info "-----Generating summary stats file for ${sampleID}-----"

  """
  python ${params.stats_agg} ${sampleID}_summary_stats.txt ${filter_stats} ${PicardMet} ${CoverageMet}
  """
}
