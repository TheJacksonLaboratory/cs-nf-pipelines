process BCF_INTERSECTVCFS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/biocontainers/bcftools:1.15--h0ea216a_2'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: confirmed_candidates_vcf

  script:

  """
  bcftools \
  isec \
  -w 1 \
  -c none \
  -n =2 \
  --threads ${task.cpus} \
  ${sampleID}_lancet.merged.vcf \
  ${sampleID}_split.vcf \
  > ${sampleID}_confirmed_candidates.vcf
"""
}