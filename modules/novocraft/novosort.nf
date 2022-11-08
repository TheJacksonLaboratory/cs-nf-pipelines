process NOVOSORT_markDuplicates {
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/novosort:lastest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'novosort' }", pattern:"*_fixed_mate_dup_marked.bam", mode:'copy'

  input:
  tuple val(sampleID), file(fixed_mate_bam)

  output:
  tuple val(sampleID), file("*_fixed_mate_dup_marked.bam"), emit: fixed_mate_dup_marked_bam

  script:
  log.info "----- Novosort Running on: ${sampleID} -----" 
  """
  novosort -markduplicates -t . -m 8G \
  ${fixed_mate_bam} > ${sampleID}_fixed_mate_dup_marked.bam
  """
