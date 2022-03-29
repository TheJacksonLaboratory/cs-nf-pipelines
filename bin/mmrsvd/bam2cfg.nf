process BAM2CFG{
  tag "$sampleID"

  cpus = 1
  time = '0:10:00'
  memory = '1.GB'
  clusterOptions = '-q batch'
  maxRetries = 1
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/bcftools:1.10.2--h4f4756c_3'

  input:
  tuple val(sampleID), file(bam)
  
  output:
  tuple val(sampleID), file(*.sv), emit: sv

  """
  bam2cfg.pl ${bam} > ${sampleID}_breakdancer.config
  breakdancer-max -r 5 -s 50 -h ${sampleID}_breakdancer.config > ${sampleID}_breakdancer.sv
  """
}
