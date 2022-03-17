process SNIFFLES {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'biocontainers/sniffles:v1.0.11ds-1-deb_cv1'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- SNIFFLES Running on: ${sampleID} -----"

  """
  sniffles -m ${bam} -v ${sampleID}_sniffles_calls.vcf
  """
}
