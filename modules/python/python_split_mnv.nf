process SPLIT_MNV {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: split_mnv_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/split_mnv.py \
  ${sampleID}_split.vcf \
  ${sampleID}_split_mnvs.vcf \
  ${tool} 
  """
}