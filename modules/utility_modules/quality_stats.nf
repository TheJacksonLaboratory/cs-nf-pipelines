process QUALITY_STATISTICS {

  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '24:00:00'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'quality_stats' }", pattern: "*fastq.gz_stat", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.fastq.gz_stat"), emit: quality_stats
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:
  log.info "----- Quality Stats Running on: ${sampleID} -----"

  if (params.read_type == "SE"){
    mode_HQ="-S -M"
    inputfq="${fq_reads[0]}"
  }
  if (params.read_type == "PE"){
    mode_HQ="-M"
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
  }

  """
  python ${projectDir}/bin/shared/filter_trim.py $mode_HQ ${params.min_pct_hq_reads}  $inputfq
  """
}
