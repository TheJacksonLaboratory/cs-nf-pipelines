process SORT {
  tag "$sampleID"

  cpus 4
  memory 20.GB
  time '20:00:00'

  container 'quay.io/jaxcompsci/samtools_with_bc:1.3.1'

  input:
  tuple val(sampleID), file(sam_file)
  val(options)

  output:
  // cannot use emit here as -n (name sort) option is incompatible with samtools index. 
  tuple val(sampleID), file("*.sorted.bam*")

  script:
  log.info "----- Samtools sort Running on: ${sampleID} -----"

  prefix = "${sampleID}.Lb"
  if (params.workflow == "chipseq"){
    output = "${prefix}.sorted.bam"
  }
  else{
    output = "${sampleID}.sorted.bam"
  }

  // check if not sorting by name
  if(options != "-n ")
  """
  samtools sort \
  ${options} \
  -@ $task.cpus \
  -O bam \
  -o ${output} \
  ${sam_file[0]}

  samtools index \
  ${output}
  """
  else
  """
  samtools sort \
  ${options} \
  -@ $task.cpus \
  -O bam \
  -o ${output} \
  ${sam_file[0]}

  """
}
