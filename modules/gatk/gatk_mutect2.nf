process GATK_MUTECT2 {
  tag "$sampleID"

  cpus = 4
  memory = 15.GB
  time {15.hour * task.attempt}
  errorStrategy 'retry' 
  maxRetries 1

  container 'broadinstitute/gatk:4.2.4.1'

  publishDir "${params.pubdir}/${ params.organize_by=='sample' ? sampleID : 'gatk' }", pattern: "*_somatic.vcf.gz", mode:'copy', enabled: params.keep_intermediate

  input:
  tuple val(sampleID), val(meta), file(normal_bam), file(normal_bai), val(normal_name), file(tumor_bam), file(tumor_bai), val(tumor_name), path(interval), val(interval_index)

  output:
  tuple val(sampleID), file("*_somatic.vcf.gz"), emit: vcf
  tuple val(sampleID), file("*_somatic.vcf.gz.tbi"), emit: tbi
  tuple val(sampleID), file("*.stats"), emit: stats

  script:
  //Estimate somatic variants using Mutect2
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus}" Mutect2 \
    -R ${params.ref_fa} \
    -I ${tumor_bam} \
    -tumor ${tumor_name} \
    -I ${normal_bam} \
    -normal ${normal_name} \
    -L ${interval} \
    --native-pair-hmm-threads 4 \
    -O ${meta.patient}_${interval_index}_somatic.vcf.gz
  """
}