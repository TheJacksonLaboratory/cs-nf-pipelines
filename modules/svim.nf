process SVIM_ALIGNMENT {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/svim:1.4.2--py_0'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- SVIM Alignment Running on: ${sampleID} -----"

  """
  svim alignment ${sampleID}_svim ${bam} ${params.fasta}
  """
}
