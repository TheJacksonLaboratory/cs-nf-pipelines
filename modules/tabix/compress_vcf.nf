process COMPRESS_VCF {
  tag tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/quay.io/biocontainers/tabix:1.11--hdfd78af_0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf.gz"), emit: compressed_vcf

  """
  bgzip \
  -c \
  ${vcf}.gz \
  > ${vcf}.gz
  """
}
