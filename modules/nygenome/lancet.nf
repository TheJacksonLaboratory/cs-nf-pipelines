process LANCET {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/lancet:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'lancet' }", pattern:".vcf", mode:'copy'

  input:
  tuple val(sampleID), file(tumor_bam)
  tuple val(sampleID), file(normal_bam)

  output:
  tuple val(sampleID), file("*_lancet.vcf"), emit: lancet_vcf

  script:
  log.info "----- Lancet Running on: ${sampleID} -----"

  """
  lancet \ 
  --tumor ${tumor_bam} \
  --normal ${normal_bam} \
  --ref ${params.ref} \
  --num-threads ${task.cpus} > ${sampleID}_lancet.vcf
  """
}
