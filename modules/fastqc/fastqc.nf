
process FASTQC {
  tag "$sampleID"

  cpus 8
  memory 4.GB
  time '10:00:00'

  container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'fastqc' }", pattern: "*_fastqc.{zip,html}", mode:'copy'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*_fastqc.{zip,html}"), emit: quality_stats


  script:

  """
    fastqc --quiet -t ${task.cpus} ${fq_reads}
  """
}
