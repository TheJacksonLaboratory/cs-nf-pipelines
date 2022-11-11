process PICARD_CLEANSAM {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*_cleaned.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_cleaned.bam"), emit: cleaned_bam

  script:

  """
  picard -Xmx24G CleanSam \
  I=${bam} \
  O=${sampleID}_cleaned.bam
  """
}
