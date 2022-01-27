/*
process SNPSIFT_DBNSFP{
container 'quay.io/biocontainers/snpsift:4.2--hdfd78af_5'

}
*/

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
