process GATK_VARIANTFILTRATION {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '03:00:00'

  container 'broadinstitute/gatk:4.2.4.1'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)
  tuple val(sampleID), file(idx)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK VariantFiltration Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  if (indel_snp == 'INDEL'){
    fs='200.0'
    output_suffix = 'INDEL_filtered.vcf'
  }
  if (indel_snp =='SNP'){
    fs ='60.0'
    output_suffix = 'SNP_filtered.vcf'
  }
  if (indel_snp == 'BOTH'){
    fs = '60.0'
    output_suffix = 'snp_indel_filtered.vcf'
  }

  """
  gatk --java-options "-Xmx${my_mem}G" VariantFiltration \
  -R ${params.ref_fa} \
  -V ${vcf} \
  -O ${sampleID}_variantfiltration_${output_suffix} \
  --cluster-window-size 10 \
  --filter-name "LowCoverage" --filter-expression "DP < 25" \
  --filter-name "VeryLowQual" --filter-expression "QUAL < 30.0" \
  --filter-expression "QUAL > 30.0 && QUAL < 50.0" --filter-name "LowQual" \
  --filter-name "LowQD" --filter-expression "QD < 1.5" \
  --filter-name "StrandBias" --filter-expression "FS > ${fs}"
  """
}