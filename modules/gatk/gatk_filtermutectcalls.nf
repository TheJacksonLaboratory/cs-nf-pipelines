process GATK_FILTERMUECTCALLS {
  tag "$sampleID"

  cpus = 1
  memory = 15.GB
  time '05:00:00'
  // time {3.hour * task.attempt}
  // errorStrategy 'retry' 
  // maxRetries 1

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_mutect2_somatic.filtered.vcf.gz", mode:'copy'
  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID+'/stats'  : 'gatk' }", pattern: "*.filteringStats.tsv", mode:'copy'

  input:
  tuple val(sampleID), file(vcf), file(tbi), file(stats)

  output:
  tuple val(sampleID), file("*_mutect2_somatic.filtered.vcf.gz"), emit: vcf
  tuple val(sampleID), file("*_mutect2_somatic.filtered.vcf.gz.tbi"), emit: tbi
  tuple val(sampleID), file("*.filteringStats.tsv"), emit: stats

  script:
  //Estimate somatic variants using Mutect2
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]
  """
  gatk --java-options "-Xmx${my_mem}G" FilterMutectCalls \
    -R ${params.ref_fa} \
    -V ${vcf} \
    --stats ${stats} \
    -O ${sampleID}_mutect2_somatic.filtered.vcf.gz
  """
}