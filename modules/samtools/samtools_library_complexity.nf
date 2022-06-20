process LIBRARY_COMPLEXITY {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(sorted_marked_bam)
  tuple val(sampleID), file(bai_file)

  output:
  tuple val(sampleID), file("tmp.ba*")

  script:
  log.info "----- Library Complexity on ${sampleID} -----"
  """
  samtools sort \
  -@ $task.cpus \
  -n \
  -O BAM \
  -o tmp.bam \
  ${sorted_marked_bam[0]}
  """
}
