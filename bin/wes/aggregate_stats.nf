process AGGREGATE_STATS_MOUSE { 
  tag "sampleID"

  cpus = 1
  time = '00:30:00'
  clusterOptions = '-q batch'

  container 'python_2.7.3.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'aggregate_stats' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(filter_stats)
  tuple val(sampleID), file(algn_met)

  output:
  tuple val(sampleID), file("*summary_stats.txt"), emit: txt

  script:
  log.info "----- Generating Summary Stats for: ${sampleID} -----"

  """
  python ${params.stats_agg} ${sampleID}_summary_stats.txt ${filter_stats} ${algn_met}
  """
}

/*
process AGGREGATE_STATS_HUMAN { D
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
*/
