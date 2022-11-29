process STRELKA2 {
  tag "$meta.patient"

  cpus = 1
  memory = 8.GB
  time = '03:00:00'

  container 'quay.io/jaxcompsci/strelka2:latest'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'strelka' }", pattern:".vcf.gz", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(meta)
  tuple val(sampleID), file(manta_vcf)

  output:
  tuple val(sampleID), file("*.indels.vcf.gz"), emit: strelka_indel_vcf
  tuple val(sampleID), file("*.snv.vcf.gz"), emit: strelka_snv_vcf

  script:
  log.info "----- Strelka2 Running on: ${sampleID} -----"

  """
  # configure strelka
  ./configureStrelkaGermlineWorkflow.py \
  --normalBam ${normal_bam} \
  --tumorBam ${tumor_bam} \
  --referenceFasta ${params.ref_fa} \
  --indelCandidates  ${projectDir}/manta/results/${sampleID}/results/variants/candidateSmallIndels.vcf.gz \
  --rundir ${sampleID}

  # execute strelka
  ${sampleID}/runWorkflow.py -m local -j ${task.cpus}
  """
}