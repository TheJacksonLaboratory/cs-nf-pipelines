process GATKv3_5_VariantRecalibrator {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'
  
  container 'broadinstitute/gatk3:3.5-0'
  
  
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.*recal"), emit: recal
  tuple val(sampleID), file("*.*trances"), emit: trances
  tuple val(sampleID), file("*.*plot.R"), emit: plot.R

  script:
  log.info "----- GATK v3.5 VariantRecalibrator Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  
  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar GenomeAnalysisTK.jar \
  -T VariantRecalibrator \
  -R ${params.ref_fa} \
  -input ${sampleID}_variants_raw.vcf \
  -resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.b37.sites.vcf \
  -resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.b37.sites.vcf \
  -resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.vcf \
  -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 dbsnp_135.b37.vcf \
  an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff \
  -mode SNP \
  -tranche 99.6
  -recalFile ${sampleID}.recal \
  -tranchesFile ${sampleID}.tranches \
  -rscriptFile ${sampleID}.plots.R
  """
}