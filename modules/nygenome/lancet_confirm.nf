process LANCET_CONFIRM {
  tag "$sampleID"

  cpus = 8
  memory = 15.GB
  time = '10:00:00'

  container 'quay.io/jaxcompsci/lancet:v1.1.0'
  // publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'lancet' }", pattern:".vcf", mode:'copy'

  input:
  tuple val(sampleID), file(bed), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name), val(chrom)

  output:
  tuple val(sampleID), file("*.vcf"), val(meta), val(normal_name), val(tumor_name), val(chrom), emit: vcf

  script:
  """
  lancet \
  --tumor ${tumor_bam} \
  --normal ${normal_bam} \
  --bed ${bed} \
  --ref ${params.ref_fa} \
  --min-k 11 \
  --low-cov 1 \
  --min-phred-fisher 5 \
  --min-strand-bias 1 \
  --min-alt-count-tumor 3 \
  --min-vaf-tumor 0.04 \
  --padding 250 \
  --window-size 2000 \
  --num-threads ${task.cpus} \
  > ${sampleID}_lancet_merged_${chrom}.vcf
  """
}
