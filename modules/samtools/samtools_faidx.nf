process SAMTOOLS_FAIDX {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/genome", mode: 'copy'

  input:
  file(fasta)

  output:
  file("*.fai")

  script:
  log.info "----- Samtools Faidx Running -----"

  """
    samtools faidx ${fasta}
  """
}
