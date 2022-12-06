process LANCET {
  tag "$meta.patient"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/lancet:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'lancet' }", pattern:".vcf", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(meta)

  output:
  tuple val(sampleID), file("*_lancet.vcf"), emit: lancet_vcf

  script:
  log.info "----- Lancet Running on: ${sampleID} -----"

  """
  lancet \ 
  --tumor ${tumor_bam} \
  --normal ${normal_bam} \
  --ref ${params.ref_fa} \
  --num-threads ${task.cpus} > ${sampleID}_lancet.vcf
  """
}
