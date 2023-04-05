process BEDTOOLS_STARTCANDIDATES {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'

  container 'quay.io/biocontainers/bedtools:2.23.0--h5b5514e_6'

  input:
  tuple val(sampleID), file(vcf), val(meta), val(chrom)

  output:
  tuple val(sampleID), file("*.vcf"), val(meta), val(chrom), emit: vcf

  script:

  """
  bedtools \
  intersect \
  -header \
  -a ${vcf} \
  -b ${params.intervalListBed} \
  -v \
  > ${sampleID}_startCand_merged_${chrom}.vcf
"""
}