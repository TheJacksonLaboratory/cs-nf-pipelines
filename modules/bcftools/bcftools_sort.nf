process BCF_SORT {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- BCFTools Sort Running on: ${sampleID} -----"

  """
  bcftools sort -o ${sampleID}_only_${indel_snp}.vcf ${vcf}
  """

}
