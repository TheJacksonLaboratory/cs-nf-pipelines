process PICARD_COLLECTRNASEQMETRICS {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '03:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.bam", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats' : 'picard' }", pattern: "*.pdf", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*metrics.txt"), emit: picard_metrics

  script:

  if (params.read_prep == "reverse_stranded") {
    strand_setting = "SECOND_READ_TRANSCRIPTION_STRAND"
  }

  if (params.read_prep == "forward_stranded") {
    strand_setting = "FIRST_READ_TRANSCRIPTION_STRAND"
  }

  if (params.read_prep == "non_stranded") {
    strand_setting = "NONE"
  }

  """
  picard CollectRnaSeqMetrics \
  I=${bam} \
  O=${sampleID}_picard_aln_metrics.txt \
  REF_FLAT=${params.ref_flat} \
  RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
  STRAND=${strand_setting} \
  CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
  """
}