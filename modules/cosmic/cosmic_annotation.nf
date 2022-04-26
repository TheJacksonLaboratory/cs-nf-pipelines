process COSMIC_ANNOTATION {
  tag "$sampleID"

  cpus 1
  memory { 10.GB * task.attempt }
  time {1.hour * task.attempt}
  clusterOptions = '-q batch'
  errorStrategy 'retry'
  maxRetries 1

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- Cosmic Annotation Running on: ${sampleID} -----"
  """
  ${projectDir}/bin/shared/Cosmic_Annotation_hg38.pl \
  -i1 ${params.cosmic} \
  -i2 ${vcf} > ${sampleID}_cosmic_annotation.vcf
  """
}
