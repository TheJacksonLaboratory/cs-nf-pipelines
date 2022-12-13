process BICSEQ2_normal {
  tag "$meta.patient"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/bicseq2:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'biqseq2' }", pattern:".txt", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)

  output:
  tuple val(sampleID), file("*.normal.params.out"), emit: normal_sample_gam

  script:

  """
  NBICseq-norm.pl \
  -l=${readLength} \
  -s=${medianInsertSize} \
  -fig=${sampleID}.GCvsRD.normal.pdf \
  -tmp=${sampleID}.tmp \
  ${configFilePathNormal} \
  ${paramsPathNormal}
  """