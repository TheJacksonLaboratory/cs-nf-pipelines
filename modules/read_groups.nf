process READ_GROUP {

  tag "sampleID"

  cpus 1
  memory 5.GB
  time '02:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*read_group.txt", mode: 'copy'

  input:
  tuple sampleID, file(read)

  output:
  tuple sampleID, file("*.txt")

  script:
  log.info "----- Read Group Information Determination Running on: ${sampleID} -----"

  """
  python ${params.read_grp_det} -p -o ${sampleID}_read_group.txt ${read[0]}
  """
  }
