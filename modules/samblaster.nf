process SAMBLASTER {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/lumpy-ref_data:0.3.1--2'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*_lumpy_aligned_bwa.bam"), emit: aligned
  tuple val(sampleID), file("*_lumpy_discordant.bam"), emit: discordant

  script:
  log.info "----- Lumpy Mapping Running on: ${sampleID} -----"

  """
  # manual read group info
  samtools view -h ${sampleID}_lumpy.bam \
  | samblaster --acceptDupMarks --excludeDups --addMateTags \
  --ignoreUnmated --maxSplitCount 2 --minNonOverlap 20 \
  | samtools view -@ ${task.cpus} -S -b - > ${sampleID}_lumpy_aligned_bwa.bam
  """
}
