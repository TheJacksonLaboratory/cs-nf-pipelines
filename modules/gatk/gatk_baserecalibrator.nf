process GATK_BASERECALIBRATOR {
  tag "$sampleID"

  cpus = 1
  memory = 35.GB
  time = '12:00:00'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'gatk' }", pattern: "*.table", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.table"), emit: table

  script:
  log.info "----- GATK BaseRecalibrator Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" BaseRecalibrator \
  -I ${bam} \
  -R ${params.ref_fa} \
  --known-sites ${params.dbSNP} \
  --known-sites ${params.gold_std_indels} \
  --known-sites ${params.phase1_1000G} \
  -O ${sampleID}_recal_data.table \
  """
}