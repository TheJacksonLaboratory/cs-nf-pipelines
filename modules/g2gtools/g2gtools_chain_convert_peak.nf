process CHAIN_CONVERT_PEAK {
  tag "$sampleID"

  cpus 1
  memory 10.GB
  time '10:00:00'

  container 'quay.io/jaxcompsci/g2gtools:0.1.31'

  input:
  tuple val(sampleID), file(bam_shifted)

  output:
  tuple val(sampleID), file("*.tmp.mm10.ba*")

  when: params.chain != null

  script:
  log.info "----- Converting Peak Coordinates to Reference on ${sampleID} -----"
  """
  g2gtools convert \
  -r -f bam -c ${params.chain} \
  -i ${bam_shifted[0]} \
  -o ${sampleID}.tmp.mm10.bam
  """
}
