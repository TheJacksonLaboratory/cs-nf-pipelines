process BCF_STARTCANDIDATES {
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

  """
  bcftools \
  intersect \
  -header \
  -a ${vcf} \
  -b ${intervalListBed} \
  -v \
  > ${sampleID}_start_candidates.vcf
"""
}