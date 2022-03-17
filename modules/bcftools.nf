process BCF_SORT {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'biocontainers/bcftools'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- BCFTools Sort Running on: ${sampleID} -----"

  """
  bcftools sort -o ${sampleID}_only_${indel_snp}.vcf ${vcf}
  """

}
process BCFTOOLS_REHEADER{
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'biocontainers/bcftools'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file(vcf), emit: vcf

  script:
  """
  bcftools reheader --samples rehead_breakdancer.txt \
  -o ${sample_name}_BreakDancerSortVCF.vcf \
  ${vcf}
  """
}
