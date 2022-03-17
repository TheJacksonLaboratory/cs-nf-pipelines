process BREAKDANCER_MAX {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/breakdancer:1.4.5--h25a10a7_7'


  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- BreakDancer Max Running on: ${sampleID} -----"

  """
  bam2cfg.pl ${bam} > ${sampleID}_breakdancer.config
  breakdancer-max -r 5 -s 50 -h ${sampleID}_breakdancer.config > ${sampleID}_breakdancer.vcf
  """
}
process FORMAT_BREAKDANCER {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'docker://python'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file('*_b2v_sorted.vcf'), emit: vcf

  script:
  log.info "----- Formatting BreakDancer on: ${sampleID} -----"

  """
  breakdancer2vcfHeader.py -i ${vcf} -o ${sampleID}_b2v.vcf
  vcfSort.sh ${sampleID}_b2v.vcf ${sampleID}_b2v_sorted.vcf
  """
}
