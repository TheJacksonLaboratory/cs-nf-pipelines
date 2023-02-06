process RENAME_VCF {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'python:3.8.10'

  input:
  tuple val(sampleID), val(meta), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: rename_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/rename_vcf.py \
  ${sampleID}_merge_prep.vcf \
  --out ${sampleID}_rename.vcf \
  ${normal}
  ${tumor}
  ${tool} 
  """
}