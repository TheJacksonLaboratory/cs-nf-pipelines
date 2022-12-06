process MANTA {
  tag "$meta.patient"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/manta:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'manta' }", pattern:".vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(meta)

  output:
  tuple val(sampleID), file("*.candidateSmallIndels.vcf.gz"), emit: manta_vcf

  script:
  log.info "----- Manta Running on: ${sampleID} -----"

  """
  # configure manta
  ./configManta.py \
  --normalBam ${normal_bam} \
  --tumorBam ${tumor_bam} \
  --referenceFasta ${params.ref_fa} \
  --callRegions ${callRegions.table} \
  --rundir ${sampleID}

  # execute manta
  ${sampleID}/runWorkflow.py -j ${task.cpus} \
  --mode local \
  --memGb ${my_mem}
  """
}