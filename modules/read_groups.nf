process READ_GROUPS {
  tag "sampleID"

  cpus 1
  memory 5.GB
  time '01:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'read_groups' }", pattern: "*read_group.txt", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  tuple val(sampleID), file("*.txt"), emit: read_groups

  script:
  log.info "----- Read Group Information Determination Running on: ${sampleID} -----"

  """
  python ${params.read_group_pyfile} -o ${sampleID}_read_group.txt ${fq_reads[0]}
  """
  }
