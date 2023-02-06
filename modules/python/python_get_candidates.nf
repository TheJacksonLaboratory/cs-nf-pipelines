process GET_CANDIDATES {
  tag "$sampleID"

  cpus 1
  memory 4.GB
  time '04:00:00'

  container 'python:3.8.10'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: rename_vcf

  script:
  """
   python \
  ${projectDir}/bin/sv/get_candidates.py \
  ${sampleID}_start_candidates.vcf\
  ${sampleID}_rename.vcf \
  ${sampleID}_candidate_merged_${chrom}.vcf
  """
}