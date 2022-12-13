process BICSEQ2_seg {
  tag "$meta.patient"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/bicseq2:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'biqseq2' }", pattern:".txt", mode:'copy'

  input:

  output:
  tuple val(sampleID), file("*.bicseq2.png"), emit: bicseq2_png
  tuple val(sampleID), file("*.bicseq2.txt"), emit: bicseq2_output

  script:

  """
  NBICseq-seg.pl \
  --control \
  --tmp ${sampleID}.tmp \
  --fig=${bicseq2PngPath} \
  --title=${sampleID} \
  --lambda=4 \
  ${segConfigFilePath} \
  ${bicseq2Path}
  """
