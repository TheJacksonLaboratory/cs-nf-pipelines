process FASTQ_SORT {

  tag "$sampleID"

  cpus 1
  memory { 50.GB * task.attempt }
  time { 2.h * task.attempt }
  errorStrategy 'retry'
  maxRetries 1

  container 'quay.io/biocontainers/fastq-tools:0.8.3--hbd632db_2'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID + '/deconvoluted_reads': 'deconvoluted_reads' }", pattern: "*.fastq", mode:'copy'


  input:
  tuple val(sampleID), file(reads)
  val(suffix)

  output:
  tuple val(sampleID), file("*sorted*{1,2}.fastq"), emit: sorted_fastq

  script:
  command_two = params.read_type == 'PE' ? "fastq-sort --id ${reads[1]} > ${sampleID}_sorted_${suffix}_2.fastq" : ''

  """
  fastq-sort --id ${reads[0]} > ${sampleID}_sorted_${suffix}_1.fastq
  ${command_two}
  """
}
