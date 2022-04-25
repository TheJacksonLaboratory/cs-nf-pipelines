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
  log.info "----- Collect RNA Sequence Metrics on: ${sampleID} -----"
  if (params.read_prep == "stranded")

    """
    picard CollectRnaSeqMetrics \
    I=${bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """

  else if (params.read_prep != "stranded")

    """
    picard CollectRnaSeqMetrics \
    I=${reordered_sorted_bam} \
    O=${sampleID}_picard_aln_metrics.txt \
    REF_FLAT=${params.ref_flat} \
    RIBOSOMAL_INTERVALS=${params.ribo_intervals} \
    STRAND=NONE \
    CHART_OUTPUT=${sampleID}_coverage_vs_transcript_plot.pdf
    """
}