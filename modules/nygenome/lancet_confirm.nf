process LANCET {
  tag "$meta.patient"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/lancet:v1.1.0'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? "$meta.patient" : 'lancet' }", pattern:".vcf", mode:'copy'

  input:
  tuple val(sampleID), file(normal_bam), file(normal_bai), val(meta)
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(meta)

  output:
  tuple val(sampleID), file("*_lancet.merged.vcf"), emit: lancet_merge_vcf

  script:

  """
  lancet \ 
  --tumor ${tumor_bam} \
  --normal ${normal_bam} \
  --bed ${chrom_bed} \
  --ref ${params.ref_fa} \
  --min-k 11 \
  --low-cov 1 \
  --min-phred-fisher 5 \
  --min-strand-bias 1 \
  --min-alt-count-tumor 3 \
  --min-vaf-tumor 0.04 \
  --padding 250 \
  --window-size 2000
  --num-threads ${task.cpus} \
  > ${sampleID}_lancet.merged.vcf
  """
}