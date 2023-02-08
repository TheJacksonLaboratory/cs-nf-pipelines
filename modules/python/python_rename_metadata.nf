process RENAME_METADATA {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: rename_metadata_vcf

  script:
  """
  python \
  ${projectDir}/bin/sv/rename_metadata.py \
  ${vcf} \
  ${vcf.baseName}_rename_metadata.vcf \
  ${tool} 
  """
}