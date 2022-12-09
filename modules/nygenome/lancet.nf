process LANCET {
  tag "$meta.patient"

  cpus = 4
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/lancet:v1.1.0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'lancet' }", pattern:"*.vcf", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name), path(bed), val(index)

  output:
  tuple val(sampleID), file("*_lancet.vcf"), emit: lancet_vcf

  script:
  """
  lancet \
  --tumor ${tumor_bam} \
  --normal ${normal_bam} \
  --ref ${params.ref_fa} \
  --bed ${bed} \
  --min-k 11 \
  --low-cov 1 \
  --min-phred-fisher 5 \
  --min-strand-bias 1 \
  --min-alt-count-tumor 3 \
  --min-vaf-tumor 0.04 \
  --num-threads ${task.cpus} > ${sampleID}_${index}_lancet.vcf
  """
}