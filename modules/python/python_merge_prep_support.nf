process MERGE_PREP_SUPPORT {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: merge_prep_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/merge_prep.py \
  --vcf ${sampleID}_rename_metadata.vcf \
  --out ${sampleID}_merge_prep.vcf \
 --tool ${tool} \
 --support
  """
}