process GATK_HAPLOTYPECALLER {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time = '10:00:00'
  clusterOptions = '-q batch'

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*.*vcf", mode:'copy'

  input:
  tuple val(sampleID), file(bam)
  tuple val(sampleID), file(bai)
  val(gvcf)

  output:
  tuple val(sampleID), file("*.*vcf"), emit: vcf
  tuple val(sampleID), file("*.idx"), emit: idx

  script:
  log.info "----- GATK Haplotype Caller Running on: ${sampleID} -----"
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  if (gvcf=='gvcf'){
    delta="-ERC GVCF"
    output_suffix='gvcf'
  }
  else{
    delta="--dbsnp ${params.dbSNP} "
    output_suffix='vcf'
  }


  """
  gatk --java-options "-Xmx${my_mem}G" HaplotypeCaller  \
  -R ${params.ref_fa} \
  -I ${bam} \
  -O ${sampleID}_variants_raw.${output_suffix} \
  -L ${params.target_gatk} \
  -stand-call-conf ${params.call_val} \
  ${params.ploidy_val} \
  ${delta} \
  """
}