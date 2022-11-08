process PICARD_FIX_MATE_INFORMATION {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'picard' }", pattern: "*fixed_mate.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*_fixed_mate.bam"), emit: fixed_mate_bam

  script:

  """
  picard -Xmx24G FixMateInformation \
  I=${bam} \
  O=${sampleID}_fixed_mate.bam \
  ADD_MATE_CIGAR=true
  """
}
