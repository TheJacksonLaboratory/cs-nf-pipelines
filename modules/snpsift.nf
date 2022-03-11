process SNPSIFT_DBNSFP{
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  // SNPEFF and SNPSIFT need updating
  container 'quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.vcf", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)
  val(indel_snp)

  output:
  tuple val(sampleID), file("*.vcf"), emit: vcf

  script:
  log.info "----- snpSift DBNSFP Running on: ${sampleID} -----"

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  if (indel_snp == 'INDEL'){
    output_suffix = 'INDEL_snpsift_dbNSFPanno.vcf'
  }
  if (indel_snp =='SNP'){
    output_suffix = 'SNP_snpsift_dbNSFPanno.vcf'
  }
  if (indel_snp == 'BOTH'){
    output_suffix = 'snp_indel_snpsift_dbNSFPanno.vcf'
  }  

  """
  java -Xmx${my_mem}G -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar \
  dbnsfp -v -db ${params.dbNSFP} -noDownload -a \
  -f SIFT_score,SIFT_pred,Polyphen2_HDIV_score,MutationAssessor_score,phyloP100way_vertebrate,1000Gp3_AF,1000Gp3_AFR_AF,1000Gp3_EUR_AF,1000Gp3_AMR_AF,1000Gp3_EAS_AF,ESP6500_AA_AF,ESP6500_EA_AF \
  ${vcf} > ${sampleID}_${output_suffix}
  """
}
process SNPSIFT_EXTRACTFIELDS {
  tag "$sampleID"

  cpus = 1
  memory = 6.GB
  time = '01:00:00'
  clusterOptions = '-q batch'

  // SNPEFF and SNPSIFT need updating
  container 'quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.txt", mode:'copy'

  input:
  tuple val(sampleID), file(vcf)

  output:
  tuple val(sampleID), file("*.txt"), emit: txt

  script:
  log.info "----- snpSift DBNSFP Running on: ${sampleID} -----"
  // add suffix for snp indel both for output name

  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  java -Xmx${my_mem}G -jar /usr/local/share/snpsift-4.3.1t-3/SnpSift.jar \
  extractFields ${vcf} \
  CHROM POS REF ALT ID FILTER QUAL FILTER AF SNPEFF_FUNCTIONAL_CLASS \
  SNPEFF_GENE_NAME SNPEFF_AMINO_ACID_CHANGE \
  SNPEFF_EFFECT SNPEFF_TRANSCRIPT_ID > ${sampleID}_snpsift_finalTable.txt
  """
}
