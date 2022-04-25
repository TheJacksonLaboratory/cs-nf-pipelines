process PICARD_COLLECTWGSMETRICS {
  tag "$sampleID"

  cpus = 1
  memory = 5.GB
  time = '08:00:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'picard' }", pattern: "*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.txt"), emit: txt

  script:
  log.info "----- Collect Alignment Sumary Metrics Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G" CollectWgsMetrics \
  --INPUT ${bam} \
  --OUTPUT ${sampleID}_CollectWgsMetrics.txt \
  --REFERENCE_SEQUENCE ${params.ref_fa} \
  --VALIDATION_STRINGENCY LENIENT
  """

}