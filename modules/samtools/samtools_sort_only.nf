process SAMTOOLS_SORT {
  tag "$sampleID"

  cpus 4
  memory 20.GB
  time '20:00:00'

  container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

  input:
  tuple val(sampleID), file(sam_file)
  val(options)
  val(suffix)

  output:
  tuple val(sampleID), file("*.sorted.*"), emit: sorted_file

  script:
  """
  samtools sort \
  ${options} \
  -@ ${task.cpus} \
  -o ${sam_file.baseName}.sorted.${suffix} \
  ${sam_file}
  """
}
