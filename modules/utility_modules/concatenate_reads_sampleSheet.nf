process CONCATENATE_READS_SAMPLESHEET {

  tag "$sampleID"

  cpus 1
  memory 15.GB
  time '03:00:00'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/concatenated_reads' : 'concatenated_reads' }", pattern: "*fastq.gz", mode:'copy'

  input:
  tuple val(sampleID), val(num_lanes), val(meta), val(read_num), path(reads)

  output:
  tuple val(sampleID), val(num_lanes), val(meta), val(read_num), path("*fastq.gz"), emit: concat_fastq

  when:
  num_lanes > 1

  script:

  """
  cat $reads > ${sampleID}_${read_num}.fastq.gz
  """
}
