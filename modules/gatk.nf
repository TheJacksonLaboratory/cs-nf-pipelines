process GATK_VARIANTANNOTATOR {
  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(sample_vcf)
  tuple val(sampleID), file(snpeff_vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- GATK VariantAnnotator Running on: ${sampleID} -----"
  """
  gatk VariantAnnotator \
  -R ${params.ref_fa} \
  -A SnpEff \
  --variant ${sample_vcf} \
  --snpEffFile ${snpeff_vcf} \
  -L ${sample_vcf} \
  -o ${sampleID}_genomeanalysistk.vcf
  """
}

process GATK_MERGEVCF {
  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(snp_vcf)
  tuple val(sampleID), file(indel_vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- GATK MergeVcfs Running on: ${sampleID} -----"
  """
  gatk MergeVcfs \
  -R ${params.ref_fa} \
  -I ${snp_vcf} \
  -I ${indel_vcf} \
  -O ${sampleID}_genomeanalysistk_combined.vcf
  """

}
// part A
process GATK_DEPTHOFCOVERAGE {

  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'
  file(params.ref_fai)

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)
  val(L)

  output:
  tuple val(sampleID), file("*_gatk_temp.txt"), emit: txt

  script:
  log.info "----- GATK Depth of Coverage Running on: ${sampleID} -----"

  """
  gatk DepthOfCoverage \
  -R ${params.ref_fa} \
  --output-format TABLE \
  -O ${sampleID}_gatk_temp.txt \
  -I ${reord_sorted_bam} \
  -L  ${L} \
  --omit-per-sample-statistics \
  --omit-interval-statistics \
  --omit-locus-table \
  """
}

process GATK_BASERECALIBRATOR {
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.table", mode:'copy'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.table"), emit: table

  script:
  log.info "----- GATK BaseRecalibrator Running on: ${sampleID} -----"

  """
  gatk BaseRecalibrator \
  -I ${bam} \
  -R ${params.ref_fa} \
  --known-sites ${params.dbSNP} \
  --known-sites ${params.gold_std_indels} \
  --known-sites ${params.phase1_1000G} \
  -O ${sampleID}_recal_data.table \
  """
}

process GATK_APPLYBQSR {
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.ba*", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(table)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  tuple val(sampleID), file("*.bai"), emit: bai

  script:
  log.info "----- GATK ApplyBQSR Running on: ${sampleID} -----"

  """
  gatk ApplyBQSR \
   -R ${params.ref_fa} \
   -I ${bam} \
   --bqsr-recal-file ${table} \
   -O ${sampleID}_realigned_BQSR.bam
  """
}


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
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK Haplotype Caller Running on: ${sampleID} -----"

  if (gvcf=='gvcf'){
    delta="-ERC GVCF"
  }
  else{
    delta="--dbsnp ${params.dbSNP} "
  }

//  --dbsnp ${params.dbSNP}

  """
  gatk HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_variants_raw.vcf \
  -L ${params.target_gatk} \
  -stand-call-conf ${params.call_val} \
  ${params.ploidy_val} \
  ${delta} \
  """
}

process GATK_SELECTVARIANTS {
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
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK Selectvariants Running on: ${sampleID} -----"

  """
  gatk SelectVariants \
  -R ${params.ref_fa} \
  -V ${vcf} \
  -select-type ${indel_snp} \
  -O ${sampleID}_selectvariants_${indel_snp}.vcf
  """
}

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
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK IndexFeatureFile Running on: ${sampleID} -----"
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
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)
  tuple val(sampleID), file(idx)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.*"), emit: vcf

  script:
  log.info "----- GATK VariantFiltration Running on: ${sampleID} -----"
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
  -O ${sampleID}_variantfiltration_${indel_snp}.vcf \
  --cluster-window-size 10 \
  --filter-name "LowCoverage" --filter-expression "DP < 25" \
  --filter-name "VeryLowQual" --filter-expression "QUAL < 30.0" \
  --filter-name "LowQD" --filter-expression "QD < 1.5" \
  --filter-name "StrandBias" --filter-expression "FS > ${fs}"
  """
}
