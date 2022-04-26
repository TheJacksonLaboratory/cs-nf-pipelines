process GATK_MERGEVCF {
  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '05:00:00'
  clusterOptions '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(snp_vcf)
  tuple val(sampleID), file(indel_vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- GATK MergeVcfs Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" MergeVcfs \
  -R ${params.ref_fa} \
  -I ${snp_vcf} \
  -I ${indel_vcf} \
  -O ${sampleID}_GATKcombined.vcf
  """
}