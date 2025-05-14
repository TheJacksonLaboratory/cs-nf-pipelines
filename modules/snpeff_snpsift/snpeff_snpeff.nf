process SNPEFF{
  tag "$sampleID"

  cpus = 1
  memory = 16.GB
  time = '06:00:00'
  errorStrategy {(task.exitStatus == 140) ? {log.info "\n\nError code: ${task.exitStatus} for task: ${task.name}. Likely caused by the task wall clock: ${task.time} or memory: ${task.memory} being exceeded.\nAttempting orderly shutdown.\nSee .command.log in: ${task.workDir} for more info.\n\n"; return 'finish'}.call() : 'finish'}

  container 'quay.io/jaxcompsci/snpeff_snpsift_5.1:v5.1d'

  publishDir "${params.pubdir}/${sampleID}", pattern:"*.*", mode:'copy', enabled: params.gen_org=='mouse' ? true : params.keep_intermediate
  publishDir "${params.pubdir}/${sampleID}", pattern:"*.*", mode:'copy', enabled: params.workflow=='amplicon_generic' ? true : params.keep_intermediate

  input:
  tuple val(sampleID),file(vcf)
  val(indel_snp)
  val(output_format)

  output:
  tuple val(sampleID),file("*.vcf"), emit:vcf
  //tuple val(sampleID),file("*.html")
  // If adding back in ^ this command should be added to the java block below
  //      -s ${sampleID}_snpeff.html \
  // tuple val(sampleID),file("*")

  script:
  
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  if (indel_snp == 'INDEL'){
    output_suffix = 'INDEL_snpeff.vcf'
  }
  if (indel_snp =='SNP'){
    output_suffix = 'SNP_snpeff.vcf'
  }
  if (indel_snp == 'BOTH'){
    output_suffix = 'SNP_INDEL_filtered_annotated_final.vcf'
  }  
  if (indel_snp == 'BOTH' && params.workflow == 'amplicon_generic' ){
    output_suffix = 'mergedCallers_filtered_annotated.vcf'
  }
  
  """
  java -Djava.io.tmpdir=./ -Xmx${my_mem}G -jar /opt/snpEff/snpEff.jar \
  ${params.gen_ver} \
  -c ${params.snpEff_config} \
  -o ${output_format} \
  -noStats \
  ${vcf} > ${sampleID}_${output_suffix}
  """
}
