process SAMTOOLS_INDEX {
  tag "sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'samtools' }", pattern:"*.ba*", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val("sampleID"), file(bam)

  output:
  tuple val("sampleID"), file("*.bai"), emit: bai

  script:
  log.info "----- Samtools Index Running on: ${sampleID} -----"

    """
    samtools index ${bam}
    """
}
