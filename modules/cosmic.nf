process COSMIC_ANNOTATION {
  tag "sampleID"

  cpus = 8
  memory = 10.GB
  time = '08:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

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
