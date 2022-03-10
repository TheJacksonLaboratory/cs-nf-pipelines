process SNPEFF{
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '06:00:00'
  clusterOptions = '-q batch'

  // SNPEFF and SNPSIFT need updating
  container 'quay.io/biocontainers/snpeff:5.1--hdfd78af_1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'snpeff' }", pattern:"*.*", mode:'copy'

  input:
  tuple val(sampleID),file(vcf)
  val(indel_snp)
  val(output_format)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf
  //tuple val(sampleID),file("*.html")
  // If adding back in ^ this command should be added to the java block below
  //          -s ${sampleID}_snpeff.html \
  // tuple val(sampleID),file("*")

  script:
  log.info "----- snpEff Running on: ${sampleID} -----"
  
  if (indel_snp == 'INDEL'){
    output_suffix = 'INDEL_snpeff.vcf'
  }
  if (indel_snp =='SNP'){
    output_suffix = 'SNP_snpeff.vcf'
  }
  if (indel_snp == 'BOTH'){
    output_suffix = 'snp_indel_snpeff.vcf'
  }  

  """
  java -Djava.io.tmpdir=$TMPDIR -Xmx8g -jar /usr/local/share/snpeff-5.1-1/snpEff.jar \
  ${params.gen_ver} \
  -c ${params.snpEff_config} \
  -o ${output_format} \
  -noStats \
  ${vcf} > ${sampleID}_${output_suffix}
  """
}

process SNPEFF_ONEPERLINE {
  tag "$sampleID"

  cpus 1
  memory 2.GB
  time '00:10:00'
  clusterOptions '-q batch'

  container 'quay.io/biocontainers/snpeff:5.1--hdfd78af_1'

  input:
  tuple val(sampleID), file(vcf)
  val(indel_snp)

  output:
  tuple	val(sampleID), file("*.vcf"), emit: vcf

  script:
  if (indel_snp == 'INDEL'){
    output_suffix = 'INDEL_snpeff.vcf'
  }
  if (indel_snp =='SNP'){
    output_suffix = 'SNP_snpeff.vcf'
  }
  if (indel_snp == 'BOTH'){
    output_suffix = 'snp_indel_snpeff.vcf'
  }
  """
  cat ${vcf} | /usr/local/share/snpeff-5.1-1/scripts/vcfEffOnePerLine.pl > ${sampleID}_oneperline_${output_suffix}
  """
}
