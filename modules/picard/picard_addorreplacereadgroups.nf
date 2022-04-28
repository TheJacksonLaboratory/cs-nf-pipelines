process PICARD_ADDORREPLACEREADGROUPS {
  tag "$sampleID"

  cpus 1
  memory 8.GB
  time '06:00:00'

  container 'quay.io/biocontainers/picard:2.26.10--hdfd78af_0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/bam' : 'picard' }", pattern: "*.bam", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(read_groups)
  tuple val(sampleID), file(bam)


  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  log.info "----- Picard Add or Replace Read Groups Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  picard AddOrReplaceReadGroups \
  INPUT=${bam} \
  OUTPUT=${sampleID}_genome_bam_with_read_groups.bam \
  SORT_ORDER=coordinate \
  \$(cat $read_groups) \
  CREATE_INDEX=true
  """
}