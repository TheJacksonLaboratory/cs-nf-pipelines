process BICSEQ2_tumor {
  tag "$meta.patient"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/bicseq2:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'biqseq2' }", pattern:".txt", mode:'copy'

  input:
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(meta)

  output:
  tuple val(sampleID), file("*.tumor.params.out"), emit: tumor_sample_gam

  script:
  log.info "----- BICSeq2-norm Running on tumor: ${sampleID} -----"

  """
  NBICseq-norm.pl \
  -l=${readLength} \
  -s=${medianInsertSize} \
  -fig=${sampleID}.GCvsRD.tumor.pdf \
  -tmp=${sampleID}.tmp \
  ${configFilePathTumor} \
  ${paramsPathTumor}
  """