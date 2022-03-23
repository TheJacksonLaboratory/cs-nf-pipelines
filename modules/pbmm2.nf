process PBMM2_INDEX {
// deprecated assumed to be precomputed
  tag "${params.fasta}"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1'

  output:
  file("*.mmi"), emit: mmi

  script:
  """
  pbmm2 index ${params.fasta} ${params.fasta.baseName}.mmi
  """
}

process PBMM2_ALIGN {
  tag "$sampleID"

  cpus = 8
  time= '72:00:00'
  memory = '250.GB'
  maxRetries = 1
  clusterOptions = '-q batch'
  errorStrategy = 'retry'

  container 'quay.io/biocontainers/pbmm2:1.3.0--h56fc30b_1'

  input:
  tuple val(sampleID), file(fq_reads)

  output:
  file "${name_string}.pbmm2.aligned.bam" into pbmm2_bam_css

  script:
  log.info "----- PBMM2 Align Running on: ${sampleID} -----"

  // MOVE THIS OUT TO THE CONFIG
  if (params.pb_mode=="ccs"){
    delta='--preset CCS'
  }
  else{
    delta='--median-filter'
  }

  """
  pbmm2 align ${params.mmi} ${fq_reads[0]} ${sampleID}_pbmm2_aligned.bam $delta --sort -j ${task.cpus}
  """
  }
