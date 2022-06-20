process CHAIN_SORT_COORDS {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/samtools-atac:1.3.1'

  input:
  tuple val(sampleID), file(bam_peak_b6_mm10)

  output:
  tuple val(sampleID), file("*.tmp1.mm10.ba*")

  when: params.chain != null

  script:
  log.info "----- Sorting bam by coordinates on ${sampleID} -----"
  """
  # sort bam by coordinates
  samtools sort \
  -@ $task.cpus -O bam \
  -o ${sampleID}.tmp1.mm10.bam \
  ${bam_peak_b6_mm10[0]}

  """
}
