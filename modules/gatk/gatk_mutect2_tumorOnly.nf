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
  tuple val(sampleID), file(tumor_bam), file(tumor_bai), val(tumor_name)

  output:
  tuple val(sampleID), file("*_somatic.vcf.gz"), file("*_somatic.vcf.gz.tbi"), file("*.stats"), emit: vcf_tbi_stats

  script:
  //Estimate somatic variants using Mutect2
  String my_mem = (task.memory-1.GB).toString()
  my_mem =  my_mem[0..-4]

  """
  gatk --java-options "-Xmx${my_mem}G -XX:ParallelGCThreads=${task.cpus}" Mutect2 \
    -R ${params.ref_fa} \
    -I ${tumor_bam} \
    -tumor ${tumor_name} \
    --germline-resource ${params.exac_ref} \
    --af-of-alleles-not-in-resource 0.0000082364 \
    --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
    --genotype-germline-sites true \
    --dont-use-soft-clipped-bases false \
    -L ${params.target_gatk} \
    --native-pair-hmm-threads 4 \
    --annotation QualByDepth \
    --annotation RMSMappingQuality \
    --annotation FisherStrand \
    --annotation MappingQualityRankSumTest \
    --annotation ReadPosRankSumTest \
    --min-base-quality-score 20 \
    -O ${sampleID}_mutect2_somatic.vcf.gz
  """
}
