process MERGE_PREP {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'quay.io/jaxcompsci/bedtools-python2:2.26.0'

  input:
  tuple val(sampleID), path(vcf), path(vci), val(normal_name), val(tumor_name)
  val(tool)

  output:
  tuple val(sampleID), path("*.vcf"), emit: merge_prep_vcf

  script:
  support_call = ${tool} == 'manta' ? '--support' : ''
  """
  python \
  ${projectDir}/bin/sv/reorder_vcf.py \
  ${vcf} \
  ${vcf.baseName}_ordered.vcf \
  ${normal_name} ${tumor_name}

  python \
  ${projectDir}/bin/sv/merge_prep.py \
  --vcf ${vcf.baseName}_ordered.vcf \
  --out ${vcf.baseName}_merge_prep.vcf \
  --tool ${tool} \
  ${support_call}
  """
}