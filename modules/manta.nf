process MANTA_CALLSV {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/manta:1.6.0--py27_0'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)
  file(fasta)

 // output:
 // tuple val(sampleID), file("*.vcf*"), emit: vcf

  script:
  log.info "Calling Manta SV"

  """
  /usr/local/bin/configManta.py --runDir mantaSVOut \
  --bam ${bam} --referenceFasta ${params.ref_fa}

  ./mantaSVOut/runWorkflow.py -m local -j ${task.cpus}

  gunzip mantaSVOut/results/variants/candidateSV.vcf.gz
  """
}