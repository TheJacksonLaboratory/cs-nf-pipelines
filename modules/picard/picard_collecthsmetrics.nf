process PICARD_COLLECTHSMETRICS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'
  clusterOptions = '-q batch'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'picard' }", pattern: "*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*Metrics.txt"), emit: hsmetrics

  script:
  log.info "----- Picard CollectHsMetrics Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  picard -Xmx${my_mem}G CollectHsMetrics \
  INPUT=${bam} \
  OUTPUT=${sampleID}_CoverageMetrics.txt \
  BAIT_INTERVALS=${params.bait_picard} \
  TARGET_INTERVALS=${params.target_picard} \
  REFERENCE_SEQUENCE=${params.ref_fa} \
  VALIDATION_STRINGENCY=SILENT
  """
}