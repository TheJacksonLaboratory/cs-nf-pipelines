process CHAIN_CONVERT_PEAK_B6 {
  tag "$sampleID"

  cpus = 1

  container 'library://taihpw/collection/g2gtools-atac:0.1.31'

  input:
  tuple val(sampleID), file(bam_shifted)

  output:
  tuple val(sampleID), file("*.tmp.mm10.ba*")

  when: params.chain != null

  script:
  log.info "----- Converting Peak Coordinates to B6 on ${sampleID} -----"
  """
  g2gtools convert \
  -r -f bam -c ${params.chain} \
  -i ${bam_shifted[0]} \
  -o ${sampleID}.tmp.mm10.bam
  """
}
