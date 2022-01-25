process SAMTOOLS_FAIDX {

  tag "sampleID"

  cpus 1
  memory 8.GB
  time '12:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'samtools' }", pattern:"*fai", mode:'copy'

  input:
  file(ref_fa)

  output:
  file("*"), emit: ref_fai

  script:
  log.info "----- Creating Reference Index -----"

    """
    samtools faidx ${ref_fa}
    """
}

process samtools index ${sampleID}_realigned_BQSR.bam
