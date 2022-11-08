process STRELKA2 {
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/manta:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'manta' }", pattern:".vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), file(tumor_bam)
  tuple val(sampleID), file(normal_bam)

  output:
  tuple val(sampleID), file("*candidateSmallIndels.vcf.gz"), emit: manta_vcf

  script:
  log.info "----- Manta Running on: ${sampleID} -----"

  """
  # configure manta
  ./configManta.py \
  --normalBam ${normal_bam} \
  --tumorBam ${tumor_bam} \
  --referenceFasta ${params.ref_fasta} \
  --rundir ${sampleID}

  # execute manta
  ${sampleID}/runWorkflow.py -j ${task.cpus}
  """