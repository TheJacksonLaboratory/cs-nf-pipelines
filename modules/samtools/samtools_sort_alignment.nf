process SORT_ALIGNMENT {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(sam_file)

  output:
  tuple val(sampleID), file("*.sorted.bam"), emit: bam
  tuple val(sampleID), file("*.sorted.bam.bai"), emit: bai

  script:
  log.info "----- Samtools sort Running on: ${sampleID} -----"
  """
  samtools sort \
  -@ $task.cpus \
  -O bam \
  -o ${sampleID}.sorted.bam \
  ${sam_file}

  samtools index \
  ${sampleID}.sorted.bam
  """

}
