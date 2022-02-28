process COSMIC_ANNOTATION {
  tag "sampleID"

  cpus = 1
  memory = 1.GB
  time = '01:00:00'
  clusterOptions = '-q batch'
// retry strategy: +=5H and +5GB

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
