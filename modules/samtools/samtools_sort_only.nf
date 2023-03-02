process SAMTOOLS_SORT {
  tag "$sampleID"

  cpus 4
  memory 20.GB
  time '20:00:00'

  container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

  input:
  tuple val(sampleID), file(sam_file)
  val(options)

  output:
  tuple val(sampleID), file("*.sorted.bam*"), emit: sorted_bam

  script:
  """
  samtools sort \
  ${options} \
  -@ ${task.cpus} \
  -O bam \
  -o ${sampleID}.sorted.bam \
  ${sam_file[0]}
  """
}
