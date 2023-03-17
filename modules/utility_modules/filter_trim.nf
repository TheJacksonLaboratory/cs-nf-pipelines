process FILTER_TRIM {

  tag "$sampleID"

  cpus 1
  memory 30.GB
  time '24:00:00'

  container 'quay.io/jaxcompsci/python-bz2file:np_2.7.18'

  publishDir "${params.pubdir}/stats", pattern: "*_stat", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*_stat"), emit: quality_stats
  tuple val(sampleID), file("*filtered_trimmed"), emit: trimmed_fastq

  script:

  if (!params.fastq2){
    inputfq="-S ${fq_reads[0]}"
  }
  else {
    inputfq="${fq_reads[0]} ${fq_reads[1]}"
  }

  """
  python ${projectDir}/bin/filter_trim.py $inputfq
  """
}