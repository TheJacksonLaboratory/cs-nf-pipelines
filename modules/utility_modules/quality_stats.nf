process QUALITY_STATISTICS {

  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '24:00:00'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'quality_stats' }", pattern: "*_stat", mode:'copy'

  input:
  tuple val(sampleID), path(fq_reads)

  output:
  tuple val(sampleID), file("*_stat"), emit: quality_stats
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:

  if (params.read_type == "SE"){
    mode_HQ="-S -M"
    inputfq="${fq_reads[0]}"
  }
  if (params.read_type == "PE"){
    mode_HQ="-M"
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
  }

  """
  python ${projectDir}/bin/shared/filter_trim.py $mode_HQ ${params.min_pct_hq_reads} -p ${params.hq_pct} $inputfq
  """
}
