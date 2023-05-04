process SAMTOOLS_FAIDX {
  tag "${fasta}"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/genome_info", mode: 'copy'

  input:
  file(fasta)

  output:
  file("*.fai")

  script:

  """
    samtools faidx ${fasta}
  """
}
