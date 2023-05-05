process LANCET {
  tag "$sampleID"

  cpus = 4
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/lancet:v1.1.0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$sampleID" : 'lancet' }", pattern:"*.vcf", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), val(meta), path(normal_bam), path(normal_bai), val(normal_name), path(tumor_bam), path(tumor_bai), val(tumor_name), path(bed), val(index)

  output:
  tuple val(sampleID), path("*_lancet.vcf"), val(meta), val(normal_name), val(tumor_name), val('lancet'), emit: vcf

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