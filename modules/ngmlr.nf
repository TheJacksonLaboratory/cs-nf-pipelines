process NGMLR_MAP{
  tag "$sampleID"

  cpus = 24
  time= '7d'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q long'
  errorStrategy = 'terminate'

  container 'quay.io/biocontainers/ngmlr:0.2.7--he513fc3_2'

  input:
  tuple val(sampleID), file(fq_reads)
  output:
  tuple val(sampleID), file "${sampleID}.sam", emit: sam

  script:
  log.info "----- NGMLR MAP Running on: ${sampleID} -----"

  """
  ngmlr -t ${task.cpus} --bam-fix -r ${params.fasta} \
  -q ${fq_reads[0]} -o ${sampleID}.sam
  """
}
