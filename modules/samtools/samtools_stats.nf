process SAMTOOLS_STATS {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/samtools' : 'samtools' }", pattern:"*.flagstat", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/samtools' : 'samtools' }", pattern:"*.idxstats", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/samtools' : 'samtools' }", pattern:"*.stats", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.flagstat")
  tuple val(sampleID), file("*.idxstats")
  tuple val(sampleID), file("*.stats")

  script:
  log.info "----- Samtools Stats Running -----"

    """
    samtools flagstat ${bam[0]} > ${bam[0]}.flagstat
    samtools idxstats ${bam[0]} > ${bam[0]}.idxstats
    samtools stats ${bam[0]} > ${bam[0]}.stats

    """
}
