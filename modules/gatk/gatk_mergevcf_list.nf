process GATK_MERGEVCF_LIST {
  tag "$sampleID"

  cpus 1
  memory 10.GB
  time '05:00:00'
  clusterOptions '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(list)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK MergeVcfs Running on: ${sampleID} -----"
  // memory needs to be set explicitly
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G" MergeVcfs \
  -R ${params.ref_fa} \
  -I ${list} \
  -O ${sampleID}_GATKcombined_raw.vcf
  """
}