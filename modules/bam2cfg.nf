process BAM2CFG {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/breakdancer:1.4.5--h25a10a7_7'
  publishDir "/home/guglib/bitbucket/", pattern: "*.config", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*.config")

  script:
  log.info "----- Summary Metrics running on ${sampleID} -----"

  """
  bam2cfg.pl ${bam} > ${sampleID}_breakdancer.config
  breakdancer-max -r ${params.breakdancer_r} -s ${params.breakdancer_s} -h ${sampleID}_breakdancer.config > ${sampleID}_breakdancer.sv
  """
}
