process QUALITY_STATISTICS {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  publishDir "${params.outdir}/quality_stats", pattern: "*fastq.gz_stat", mode: 'copy'

  input:
  tuple val(sampleID), file(reads)

  output:
  tuple val(sampleID), file("*.fastq.gz_stat"), emit: quality_stats
  tuple val(sampleID), file("${sampleID}_R{1,2}*filtered_trimmed"), emit: trimmed_fastq

  script:
  log.info "----- Quality Stats Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    mode_HQ="-S -M"
    inputfq="${reads}"
  }
  if (params.read_type == "PE"){
    mode_HQ="-M"
    inputfq="${reads[0]} ${reads[1]}"
  }

  """
  python ${params.filter_trim} $mode_HQ ${params.min_pct_hq_reads}  $inputfq
  """
}
