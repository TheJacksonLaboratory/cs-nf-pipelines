
process SNPSIFT_DBNSFP{
  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  container 'gatk-4.1.6.0_samtools-1.3.1_snpEff_4.3_vcftools_bcftools.sif'

  input:
  tuple val(sampleID), file(vcf)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf
  
  script:
  """
  java -jar /snpEff_v4_3/snpEff/SnpSift.jar \
  dbnsfp -v -db ${params.dbNSFP} -noDownload -a \
  -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF \
  ${vcf} > ${sampleID}_${indel_snp}.vcf
  """
}


process SNPSIFT_EXTRACTFIELDS {

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'
  
  container 'gatk-3.6_snpeff-3.6c_samtools-1.3.1_bcftools-1.11.sif'
//container 'quay.io/biocontainers/snpsift:4.2--hdfd78af_5'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.txt"), emit: txt  

  script:
  """
  java -jar /snpEff/SnpSift.jar \
  extractFields ${vcf} \
  CHROM POS REF ALT ID FILTER QUAL FILTER AF SNPEFF_FUNCTIONAL_CLASS \
  SNPEFF_GENE_NAME SNPEFF_AMINO_ACID_CHANGE \
  SNPEFF_EFFECT SNPEFF_TRANSCRIPT_ID > ${sampleID}_snpsift.txt
  """

}
