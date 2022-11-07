process GATK_FIX_MATE_INFORMATION {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*fixed_mate.bam", mode:'copy'

  input:
  tuple val(sampleID), file(marked_bam)

  output:
  tuple val(sampleID), file("*_fixed_mate.bam"), emit: fixed_mate_bam

  script:
  log.info "----- GATK FixMateInformation Running on: ${sampleID} -----"

  """
  picard -Xmx24G FixMateInformation \
  I=${marked_bam} \
  O=${sampleID}_fixed_mate.bam \
  ADD_MATE_CIGAR=true
  """
}
