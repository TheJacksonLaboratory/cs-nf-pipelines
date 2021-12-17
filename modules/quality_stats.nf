process QUALITY_STATISTICS {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*fastq.gz_stat", mode: 'copy'

  input:
  tuple val(sampleID), file(read1), file(read2)

  output:
  file "*.fastq.gz_stat"
  tuple sampleID, file("*.fastq.gz_stat")
  tuple sampleID, file("${sampleID}_R{1,2}*filtered_trimmed")

  script:
  log.info "----- Quality Stats Running on: ${sampleID} -----"

  if (params.reads == "SE"){
    mode_HQ="-S -M"
    inputfq="${read1}"
  }
  if (params.reads == "PE"){
    mode_HQ="-M"
    inputfq="${read1} ${read2}"
  }

  """
  python ${params.filter_trim} $mode_HQ ${params.min_pct_hq_reads}  $inputfq
  """
}
