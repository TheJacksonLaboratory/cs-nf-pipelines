process SAMTOOLS_INDEX {
  tag "$sampleID"

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
process SAMTOOLS_STATS {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'samtools' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val("sampleID"), file(bam)

  output:
  tuple val("sampleID"), file("*.txt"), emit: txt

  script:
  log.info "----- Samtools Stats Running on: ${sampleID} -----"

    """
    samtools stats ${bam} |grep "^IS" |awk '{a = a + \$2*\$3; b = b + \$3}END{print int(a/b)}' > ${sampleID}_insert_size.txt
    """
}

process SAMTOOLS_VIEW {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'samtools' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val("sampleID"), file(bam)

  output:
  tuple val("sampleID"), file("*.bam"), emit: bam

  script:
  log.info "----- Samtools Stats Running on: ${sampleID} -----"

  """
  # Extract the discordant pairedend alignments
  samtools view -@ ${task.cpus} -b ${params.samtools_flag} ${bam} > ${sampleID}_discordant.bam
  """
}

