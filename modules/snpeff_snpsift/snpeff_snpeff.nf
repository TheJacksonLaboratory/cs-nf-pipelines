process SNPEFF{
  tag "$sampleID"

  cpus = 1
  memory = 8.GB
  time = '06:00:00'

  // SNPEFF and SNPSIFT need updating
  container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1'

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
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

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
  java -Djava.io.tmpdir=$TMPDIR -Xmx${my_mem}G -jar /opt/snpEff/snpEff.jar \
  ${params.gen_ver} \
  -c ${params.snpEff_config} \
  -o ${output_format} \
  -noStats \
  ${vcf} > ${sampleID}_${output_suffix}
  """
}