process READ_GROUPS {
  // ** why is only one read used here?
  tag "sampleID"

  cpus 1
  memory 5.GB
  time '02:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  publishDir "${outdir}/read_groups", pattern: "*read_group.txt", mode: 'copy'

  input:
  tuple val(sampleID), file(read)

  output:
  tuple val(sampleID), file("*.txt"), emit: read_group

  script:
  log.info "----- Read Group Information Determination Running on: ${sampleID} -----"

  """
  python ${params.read_grp_det} -p -o ${sampleID}_read_group.txt ${read[0]}
  """
  }
