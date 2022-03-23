process PBSV_DISCOVERY{
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/pbsv-td_refs:2.3.0--0'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.svsig.gz"), emit: svsig

  script:
  log.info "----- PBSV Discover Running on: ${sampleID} -----"

  if (params.tandem==true){
    delta="--tandem-repeats ${params.defultname}"
  }
  else{
    delta=""
  }

  """
  pbsv discover $delta ${bam} ${sampleID}.svsig.gz
  """
}

process PBSV_CALL{
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/jaxcompsci/pbsv-td_refs:2.3.0--0'

  input:
  tuple val(sampleID), file(svsig)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- PBSV Call Running on: ${sampleID} -----"

  if (params.pb_mode=="ccs"){
    delta='--ccs'
  }
  else{
    delta=''
  }

  """
  pbsv call $delta ${params.fasta} ${svsig} pbsv_calls.vcf
  """
}
