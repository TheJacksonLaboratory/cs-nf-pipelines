process BISMARK_DEDUPLICATION {
  tag "$sampleID"

  cpus 8
  memory {60.GB * task.attempt}
  time {30.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'quay.io/biocontainers/bismark:0.23.1--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/alignment' : 'bismark_align' }", pattern: "*.bam", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'bismark_align' }", pattern: "*txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.bam"), emit: dedup_bam
  tuple val(sampleID), file("*report.txt"), emit: dedup_report

  script:
  log.info "----- Bismark Deduplication Running on: ${sampleID} -----"

  fq_type = params.read_type == 'PE' ? '-p' : '-s'
  
  """
  deduplicate_bismark $fq_type --bam $bam
  """
}