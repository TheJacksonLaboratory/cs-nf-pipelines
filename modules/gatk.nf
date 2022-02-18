process GATK_PRINTREADS{
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(table)

  output:
  tuple val(sampleID), file("*.bam"), emit: bam
  
  script:
  log.info "----- GATK PrintReads Running on: ${sampleID} -----"
// ran into an issue with using the old -BQSR argument and using our table from baserecal
  """
  gatk PrintReads \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_realigned_BQSR.bam \
  """
}

process GATK_REALIGNERTARGETCREATOR {
  tag "sampleID"

  cpus = 12
  memory = 35.GB
  time = '72:00:00'
  clusterOptions = '-q batch'

  // update container, cannot find this command in gatk 4
  // container 'broadinstitute/gatk:4.2.4.1'

  container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'

  input:
  tuple val(sampleID), file(bam)

  output:
  tuple val(sampleID), file("*.intervals"), emit: intervals

  script:
  log.info "----- GATK RealignerTargetCreator Running on: ${sampleID} -----"

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx24g -jar /usr/GenomeAnalysisTK.jar \
  -I ${bam} \
  -R ${params.ref_fa} \
  -T RealignerTargetCreator \
  -o ${sampleID}.aligner.intervals \
  -nt 12 \
  --disable_auto_index_creation_and_locking_when_reading_rods
  """

}

process GATK_INDELREALIGNER{
tag "sampleID"

cpus = 12
memory = 35.GB
time = '72:00:00'
clusterOptions = '-q batch'

// update container, cannot find this command in gatk 4
// container 'broadinstitute/gatk:4.2.4.1'

container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'

input:
tuple val(sampleID), file(bam)
tuple val(sampleID), file(intervals)

output:
tuple val(sampleID), file("*.bam"), emit: bam
tuple val(sampleID), file("*.bai"), emit: bai

script:
log.info "----- GATK IndelRealigner Running on: ${sampleID} -----"

"""
java -Djava.io.tmpdir=$TMPDIR -Xmx24g -jar /usr/GenomeAnalysisTK.jar \
-I ${bam} \
-R ${params.ref_fa} \
-T IndelRealigner \
-targetIntervals ${intervals} \
-o ${sampleID}_realigned.bam \
--disable_auto_index_creation_and_locking_when_reading_rods
"""
}

process GATK_VARIANTANNOTATOR {
  tag "sampleID"

  cpus 1
  memory 15.GB
  time '24:00:00'
  clusterOptions '-q batch'

 // container 'broadinstitute/gatk:4.2.4.1'
  container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'

  input:
  tuple val(sampleID), file(sample_vcf)
  tuple val(sampleID), file(snpeff_vcf)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- GATK VariantAnnotator Running on: ${sampleID} -----"
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
  -T VariantAnnotator \
  -R ${params.ref_fa} \
  -A SnpEff \
  --variant ${sample_vcf} \
  --snpEffFile ${snpeff_vcf} \
  -L ${sample_vcf} \
  -o ${sampleID}_GATKannotated.vcf
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
  -O ${sampleID}_GATKcombined.vcf
  """

}
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
  -I ${bam} \
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
    output_suffix='gvcf'
  }
  else{
    delta="--dbsnp ${params.dbSNP} "
    output_suffix='vcf'
  }

//  --dbsnp ${params.dbSNP}

  """
  gatk HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_variants_raw.${output_suffix} \
  -L ${params.target_gatk} \
  -stand-call-conf ${params.call_val} \
  ${params.ploidy_val} \
  ${delta} \
  """
}

process GATK_HAPLOTYPECALLER_WGS {
  tag "sampleID"

  cpus = 8
  memory = 15.GB
  time = '23:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:

  log.info "----- GATK Haplotype Caller Running on: ${sampleID} -----"

  """
  gatk HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_HaplotypeCaller.vcf \
  -stand-call-conf ${params.call_val}
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
  tuple val(sampleID), file("*.vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK VariantFiltration Running on: ${sampleID} -----"
  if (indel_snp == 'INDEL'){
    fs='200.0'
    output_suffix = 'INDEL_filtered.vcf'
  }
  if (indel_snp =='SNP'){
    fs ='60.0'
    output_suffix = 'SNP_filtered.vcf'
  }
  if (indel_snp == 'MOUSE'){
    // mouse will be indel but fs needs to be same as snp (not sure why)
    fs = '60.0'
    output_suffix = 'snp_indel_filtered.vcf'
  }

  """
  gatk VariantFiltration \
  -R ${params.ref_fa} \
  -V ${vcf} \
  -O ${sampleID}_variantfiltration_${output_suffix} \
  --cluster-window-size 10 \
  --filter-name "LowCoverage" --filter-expression "DP < 25" \
  --filter-name "VeryLowQual" --filter-expression "QUAL < 30.0" \
  --filter-name "LowQD" --filter-expression "QD < 1.5" \
  --filter-name "StrandBias" --filter-expression "FS > ${fs}"
  """
}
