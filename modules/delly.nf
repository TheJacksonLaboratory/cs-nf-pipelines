process DELLY_CALL {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/delly-ref_data:0.8.3--hf3ca161_1'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

//  output:
//  tuple val(sampleID), file()

  script:
  log.info "----- Delly Running on: ${sampleID} -----"

  """
  delly call \
  -q 40 \
  -x ${params.exclude_regions} \
  -s 500 \
  -o ${sampleID}_delly.bcf \
  -g ${params.ref_fa} ${bam}
  """
}
