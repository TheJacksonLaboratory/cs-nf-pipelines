// part A
process GATK_STATS_A {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'
  file(params.ref_fai)

  input:
  tuple val(sampleID), file(reord_sorted_bam)
  tuple val(sampleID), file(reord_sorted_bai)

  output:
  tuple val(sampleID), file("*gatk_temp3*"), emit: gatk_3
  tuple val(sampleID), file("*gatk_temp6*"), emit: gatk_6

  when:
  params.gen_org == "human"

  script:
  log.info "----- Human GATK Coverage Stats, Part 1 Running on: ${sampleID} -----"
  """
  gatk DepthOfCoverage \
  -R ${params.ref_fa} \
  --output-format TABLE \
  -O ${sampleID}_gatk_temp1.txt \
  -I ${reord_sorted_bam} \
  -L  ${params.probes} \
  --omit-per-sample-statistics \
  --omit-interval-statistics \
  --omit-locus-table \

  gatk DepthOfCoverage \
  -R ${params.ref_fa} \
  --output-format TABLE \
  -O ${sampleID}_gatk_temp4.txt \
  -I ${reord_sorted_bam} \
  -L ${params.ctp_genes} \
  --omit-per-sample-statistics \
  --omit-interval-statistics \
  --omit-locus-table \

  chmod +x ${params.gatk_form}

  ${params.gatk_form} ${sampleID}_gatk_temp1.txt ${sampleID}_gatk_temp2.txt ${sampleID}_gatk_temp3.txt ${params.probes}

  ${params.gatk_form} ${sampleID}_gatk_temp4.txt ${sampleID}_gatk_temp5.txt ${sampleID}_gatk_temp6.txt ${params.ctp_genes}
  """
}
// part B
process GATK_STATS_B {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'python_2.7.sif'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*", mode:'copy'

  input:
  tuple val(sampleID), file(gatk_3)
  tuple val(sampleID), file(gatk_6)

  output:
  file "*CCP_interval_avg_median_coverage.bed"
  file "*exome_interval_avg_median_coverage.bed"
  tuple val(sampleID), file("*CP_interval_avg_median_coverage.bed")

  when:
  params.gen_org == "human"

  script:
  log.info "----- Human GATK Coverage Stats, Part 2 Running on: ${sampleID} -----"

  """
  python ${params.cov_calc} ${gatk_3} ${sampleID}_exome_interval_avg_median_coverage.bed

  python ${params.cov_calc} ${gatk_6} ${sampleID}_CCP_interval_avg_median_coverage.bed

  """
}
/*
process GATK_BASERECALIBRATOR {
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*realigned_BQSR*", mode: 'copy'

  input:
  tuple sampleID, file(dedup_bam)

  output:
  tuple sampleID, file("*.table") emit: recal_data_table

  script:
  log.info "----- GATK BaseRecalibrator Running on: ${sampleID} -----"

  """
  gatk BaseRecalibrator \
  -I ${dedup_bam} \
  -R ${params.ref_fa} \
  --known-sites ${params.dbSNP} \
  --known-sites ${params.gold_std_indels} \
  --known-sites ${params.phase1_1000G} \
  -O recal_data.table \
  """
}
process GATK_APPLYBQSR {
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  publishDir "${sample_tmpdir}_tmp", pattern: "*realigned_BQSR*", mode: 'copy'

  input:
  tuple sampleID, file(dedup_bam)

  output:
  tuple sampleID, file("*realigned_BQSR.bam")
  tuple sampleID, file("*realigned_BQSR*bai")

  script:
  log.info "----- GATK ApplyBQSR Running on: ${sampleID} -----"

  """
  gatk --java-options "-Xmx24g" ApplyBQSR \
   -R ${params.ref_fa} \
   -I ${dedup_bam} \
   --bqsr-recal-file recal_data.table \
   -O ${sampleID}_realigned_BQSR.bam
  """
}
*/

process GATK_HAPLOTYPECALLER {
  tag "sampleID"

  cpus = 8
  memory = 15.GB
  time = '23:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)
  val(gvcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- GATK Haplotype Caller Running on: ${sampleID} -----"

  if (gvcf=='gvcf'){
    delta='-ERC GVCF'
  }
  else{
    delta=''
  }

//  --dbsnp ${params.dbSNP}

  """
  gatk HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_variants_raw.vcf \
  -ERC GVCF \
  """
}
/*
process GATK_SELECTVARIANTS{
  tag "sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'


  input:
  tuple sampleID, file(merged_raw_var) from raw_variants

  output:
  tuple sampleID, file("*only_snps_filtered.vcf"), file("*only_indels_filtered.vcf") into filt_var, dummy_filt_var

  input:
  tuple val(sampleID), file(merged_raw_var)
  val(indel_snp)

   output:
   tuple val(sampleID), file(*.vcf)

  """
  gatk SelectVariants \
  -R ${params.ref_fa} \
  -V ${merged_raw_var} \
  -select-type ${indel_snp} \
  -O ${sampleID}_${indel_snp}.vcf
  """
}*/

process GATK_INDEXFEATUREFILE {
  tag "sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.idx", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)
 
  output:
  tuple val(sampleID), file("*.idx"), emit: vcf_index
 
  script:
  
  """
  gatk IndexFeatureFile \
  -I ${vcf}
  """
}


process GATK_VARIANTFILTRATION {
  tag "sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'
  
  input:
  tuple val(sampleID), file(vcf)
  tuple val(sampleID), file(idx)
  val(indel_snp)
  
  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  if (indel_snp == 'INDEL'){
    fs='200.0'
  }
  if (indel_snp =='SNP'){
    fs ='60.0'
  }
  if (params.gen_org == 'mouse'){
    // mouse will be indel but fs needs to be same as snp (not sure why)
    fs = '60.0'
  }

  """
  gatk VariantFiltration \
  -R ${params.ref_fa} \
  -V ${vcf} \
  -O ${sampleID}_${indel_snp}.vcf \
  --cluster-window-size 10 \
  --filter-expression "DP < 25" --filter-name "LowCoverage" \
  --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
  --filter-expression "QUAL > 30.0 && QUAL < 50.0" --filter-name "LowQual" \
  --filter-expression "QD < 1.5" --filter-name "LowQD" \
  --filter-expression "FS > ${fs}" --filter-name "StrandBias"
  """
}
