process FILTER_REMOVE_MULTI_SIEVE {
  tag "$sampleID"

  cpus 8
  memory 4.GB
  time '04:00:00'

  container 'quay.io/biocontainers/deeptools:3.3.2--py_1'

  input:
  tuple val(sampleID), file(shift_bams)

  output:
  tuple val(sampleID), file("*.shift.tmp.ba*")

  script:
  log.info "----- Running deeptools alignmentSieve on ${sampleID} -----"
  """
  alignmentSieve \
  --numberOfProcessors $task.cpus \
  --ATACshift \
  --bam ${shift_bams[0]} \
  -o ${sampleID}.shift.tmp.bam 
  """

}
