process COMPRESS_SOMATIC_VCFs {
  tag tag "$meta.patient"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/quay.io/biocontainers/tabix:1.11--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'compressed_somatic_vcf' }", pattern:".vcf.gz", mode:'copy'

  input:
  # this is a place holder
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)

  output:
  # also a place holder
  tuple val(sampleID), file("*_lancet.vcf"), emit: lancet_vcf

  """
  bgzip \
  -c \
  ${vcf} \
  > ${vcf}.gz
  """
}
