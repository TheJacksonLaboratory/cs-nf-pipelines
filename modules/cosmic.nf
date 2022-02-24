process COSMIC_ANNOTATION {
  tag "sampleID"

  cpus = 8
  memory = 10.GB
  time = '08:00:00'
  clusterOptions = '-q batch'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- Cosmic Annotation Running on: ${sampleID} -----"
  """
  ${params.cosmic_annot} \
  -i1 ${params.cosmic} \
  -i2 ${vcf} > ${sampleID}_cosmic_annotation.vcf
  """
}
