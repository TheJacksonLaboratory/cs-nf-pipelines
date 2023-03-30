process SAMTOOLS_FAIDX {
  tag "$sampleID"

  cpus 1
  memory 2.GB
  time '00:30:00'

  container 'quay.io/biocontainers/samtools:1.14--hb421002_0'

  publishDir "${params.pubdir}/g2gtools", pattern:"*.fai", mode:'copy', enabled: params.workflow == 'generate_pseudoreference' ? true : false

  input:
  tuple val(sampleID), path(fasta)

  output:
  tuple val(sampleID), path("*.fai"), emit: fai

  script:

    """
    samtools faidx ${fasta}
    """
}
