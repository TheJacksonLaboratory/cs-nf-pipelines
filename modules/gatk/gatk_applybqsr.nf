process GATK_APPLYBQSR {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'gatk' }", pattern: "*.bam", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(table)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  log.info "----- GATK ApplyBQSR Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" ApplyBQSR \
  -R ${params.ref_fa} \
  -I ${bam} \
  --bqsr-recal-file ${table} \
  -O ${sampleID}_realigned_BQSR.bam
  """
}