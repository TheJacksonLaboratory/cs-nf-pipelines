process PICARD_SORTSAM {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*_sortsam.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(sam)

  output:
  tuple val(sampleID), file("*_sortsam.bam"), emit: bam
  tuple val(sampleID), file("*_sortsam.bai"), emit: bai

  script:
  log.info "----- Picard SortSam Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  picard -Xmx${my_mem}G SortSam \
  SO=coordinate \
  INPUT=${sam} \
  OUTPUT=${sampleID}_sortsam.bam  \
  VALIDATION_STRINGENCY=SILENT \
  CREATE_INDEX=true
  """
}